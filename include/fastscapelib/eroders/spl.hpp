/**
 * Functions to compute bedrock channel erosion.
 */
#ifndef FASTSCAPELIB_ERODERS_SPL_H
#define FASTSCAPELIB_ERODERS_SPL_H

#include <cassert>
#include <cmath>
#include <limits>
#include <type_traits>

#include "xtensor/xtensor.hpp"
#include "xtensor/xmanipulation.hpp"

#include "fastscapelib/utils/utils.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"


namespace fastscapelib
{

    template <class FG, class S = typename FG::xt_selector>
    class spl_eroder
    {
    public:
        using flow_graph_type = FG;
        using xt_selector = S;
        using grid_type = typename flow_graph_type::grid_type;

        using size_type = typename grid_type::size_type;
        using shape_type = typename grid_type::shape_type;
        using data_type = typename grid_type::grid_data_type;
        using data_array_type = xt_array_t<xt_selector, data_type>;
        using elevation_type = data_array_type;

        template <class K>
        spl_eroder(FG& flow_graph, K&& k_coef, double m_exp, double n_exp, double tolerance)
            : m_flow_graph(flow_graph)
            , m_shape(flow_graph.grid().shape())
            , m_tolerance(tolerance)
        {
            set_k_coef(k_coef);
            set_m_exp(m_exp);
            set_n_exp(n_exp);

            m_erosion = xt::zeros<data_type>(m_shape);
        };

        const data_array_type& k_coef()
        {
            return m_k_coef;
        };

        template <class T>
        void set_k_coef(T value, typename std::enable_if_t<std::is_floating_point<T>::value>* = 0)
        {
            m_k_coef = xt::broadcast(std::forward<T>(value), { m_flow_graph.size() });
        };

        template <class K>
        void set_k_coef(K&& value, typename std::enable_if_t<xt::is_xexpression<K>::value>* = 0)
        {
            assert(value.shape() == m_shape);
            m_k_coef = xt::flatten(value);
        };

        double m_exp()
        {
            return m_m_exp;
        };

        void set_m_exp(double value)
        {
            // TODO: validate value
            m_m_exp = value;
        };

        double n_exp()
        {
            return m_n_exp;
        };

        void set_n_exp(double value)
        {
            // TODO: validate value
            m_n_exp = value;
        };

        double tolerance()
        {
            return m_tolerance;
        }

        size_type n_corr()
        {
            return m_n_corr;
        };

        const data_array_type& erode(const data_array_type& elevation,
                                     const data_array_type& drainage_area,
                                     double dt);

    private:
        flow_graph_type& m_flow_graph;
        shape_type m_shape;
        data_array_type m_k_coef;
        data_array_type m_erosion;
        double m_m_exp;
        double m_n_exp;
        double m_tolerance;
        size_type m_n_corr;
    };


    template <class FG, class S>
    auto spl_eroder<FG, S>::erode(const data_array_type& elevation,
                                  const data_array_type& drainage_area,
                                  double dt) -> const data_array_type&
    {
        auto erosion_flat = xt::flatten(m_erosion);
        const auto elevation_flat = xt::flatten(elevation);
        const auto drainage_area_flat = xt::flatten(drainage_area);

        auto& flow_graph_impl = m_flow_graph.impl();

        const auto& receivers = flow_graph_impl.receivers();
        const auto& receivers_count = flow_graph_impl.receivers_count();
        const auto& receivers_distance = flow_graph_impl.receivers_distance();
        const auto& receivers_weight = flow_graph_impl.receivers_weight();
        const auto& dfs_stack = flow_graph_impl.dfs_stack();

        m_n_corr = 0;

        for (const auto& istack : dfs_stack)
        {
            auto r_count = receivers_count[istack];

            for (size_type r = 0; r < r_count; ++r)
            {
                size_type irec = receivers(istack, r);

                if (irec == istack)
                {
                    // at basin outlet or pit
                    erosion_flat(istack) = 0.;
                    continue;
                }

                data_type istack_elevation = elevation_flat(istack);  // at time t
                data_type irec_elevation
                    = elevation_flat(irec) - erosion_flat(irec);  // at time t+dt

                if (irec_elevation >= istack_elevation)
                {
                    // may happen if flow is routed outside of a depression / flat area
                    erosion_flat(istack) = 0.;
                    continue;
                }

                auto factor
                    = (m_k_coef(istack) * dt * std::pow(drainage_area_flat(istack), m_m_exp));

                data_type delta_0 = istack_elevation - irec_elevation;
                data_type delta_k;

                if (m_n_exp == 1)
                {
                    // fast path for n_exp = 1 (common use case)
                    factor /= receivers_distance(istack, r);
                    delta_k = delta_0 / (1. + factor);
                }

                else
                {
                    // 1st order Newton-Raphson iterations (k)
                    factor /= std::pow(receivers_distance(istack, r), m_n_exp);
                    delta_k = delta_0;

                    while (true)
                    {
                        auto factor_delta_exp = factor * std::pow(delta_k, m_n_exp);
                        auto func = delta_k + factor_delta_exp - delta_0;

                        if (func <= m_tolerance)
                        {
                            break;
                        }

                        auto func_deriv = 1 + m_n_exp * factor_delta_exp / delta_k;
                        delta_k -= func / func_deriv;

                        if (delta_k <= 0)
                        {
                            break;
                        }
                    }
                }

                if (delta_k <= 0)
                {
                    // prevent the creation of new depressions / flat channels
                    // by arbitrarily limiting erosion
                    m_n_corr++;
                    delta_k = std::numeric_limits<data_type>::min();
                }

                if (r_count == 1)
                {
                    erosion_flat(istack) = (delta_0 - delta_k);
                }
                else
                {
                    erosion_flat(istack) += (delta_0 - delta_k) * receivers_weight(istack, r);
                }
            }
        }

        return m_erosion;
    };


    template <class FG, class K>
    spl_eroder<FG> make_spl_eroder(
        FG& graph, K&& k_coef, double m_exp, double n_exp, double tolerance)
    {
        return spl_eroder<FG>(graph, k_coef, m_exp, n_exp, tolerance);
    }


    namespace detail
    {

        template <class T, class S>
        auto k_coef_as_array(T k_coef,
                             S&& shape,
                             typename std::enable_if_t<std::is_floating_point<T>::value>* = 0)
        {
            return xt::broadcast(std::forward<T>(k_coef), shape);
        }


        template <class K, class S>
        auto k_coef_as_array(K&& k_coef,
                             S&& shape,
                             typename std::enable_if_t<xt::is_xexpression<K>::value>* = 0)
        {
            auto k_coef_arr = xt::flatten(k_coef);

            assert(k_coef_arr.shape() == shape);
            (void) shape;  // TODO: still unused parameter warning despite assert?

            return k_coef_arr;
        }


        /**
         * erode_stream_power implementation.
         *
         * The implementation here slightly differs from the one described in
         * Braun & Willet (2013). The problem is reformulated so that the
         * Newton-Raphson method is applied on the difference of elevation
         * between a node and its receiver, rather than on the node's
         * elevation itself. This allows saving some operations.
         */
        template <class Er, class El, class Dr, class FG, class K>
        auto erode_stream_power_impl(Er&& erosion,
                                     El&& elevation,
                                     Dr&& drainage_area,
                                     FG& flow_graph,
                                     K&& k_coef,
                                     double m_exp,
                                     double n_exp,
                                     double dt,
                                     double tolerance) -> typename FG::size_type
        {
            using T = std::common_type_t<typename std::decay_t<Er>::value_type,
                                         typename std::decay_t<El>::value_type>;
            using size_type = typename FG::size_type;

            auto erosion_flat = xt::flatten(erosion);
            const auto elevation_flat = xt::flatten(elevation);
            const auto drainage_area_flat = xt::flatten(drainage_area);
            const auto k_coef_arr = k_coef_as_array(k_coef, elevation_flat.shape());

            auto& flow_graph_impl = flow_graph.impl();

            const auto& receivers = flow_graph_impl.receivers();
            const auto& receivers_count = flow_graph.impl().receivers_count();
            const auto& receivers_distance = flow_graph.impl().receivers_distance();
            const auto& receivers_weight = flow_graph.impl().receivers_weight();
            const auto& dfs_stack = flow_graph.impl().dfs_stack();

            size_type n_corr = 0;

            for (const auto& istack : dfs_stack)
            {
                auto r_count = receivers_count[istack];

                for (size_type r = 0; r < r_count; ++r)
                {
                    size_type irec = receivers(istack, r);

                    if (irec == istack)
                    {
                        // at basin outlet or pit
                        erosion_flat(istack) = 0.;
                        continue;
                    }

                    T istack_elevation = elevation_flat(istack);                   // at time t
                    T irec_elevation = elevation_flat(irec) - erosion_flat(irec);  // at time t+dt

                    if (irec_elevation >= istack_elevation)
                    {
                        // may happen if flow is routed outside of a depression / flat area
                        erosion_flat(istack) = 0.;
                        continue;
                    }

                    auto factor
                        = (k_coef_arr(istack) * dt * std::pow(drainage_area_flat(istack), m_exp));

                    T delta_0 = istack_elevation - irec_elevation;
                    T delta_k;

                    if (n_exp == 1)
                    {
                        // fast path for n_exp = 1 (common use case)
                        factor /= receivers_distance(istack, r);
                        delta_k = delta_0 / (1. + factor);
                    }

                    else
                    {
                        // 1st order Newton-Raphson iterations (k)
                        factor /= std::pow(receivers_distance(istack, r), n_exp);
                        delta_k = delta_0;

                        while (true)
                        {
                            auto factor_delta_exp = factor * std::pow(delta_k, n_exp);
                            auto func = delta_k + factor_delta_exp - delta_0;

                            if (func <= tolerance)
                            {
                                break;
                            }

                            auto func_deriv = 1 + n_exp * factor_delta_exp / delta_k;
                            delta_k -= func / func_deriv;

                            if (delta_k <= 0)
                            {
                                break;
                            }
                        }
                    }

                    if (delta_k <= 0)
                    {
                        // prevent the creation of new depressions / flat channels
                        // by arbitrarily limiting erosion
                        n_corr++;
                        delta_k = std::numeric_limits<T>::min();
                    }

                    if (r_count == 1)
                    {
                        erosion_flat(istack) = (delta_0 - delta_k);
                    }
                    else
                    {
                        erosion_flat(istack) += (delta_0 - delta_k) * receivers_weight(istack, r);
                    }
                }
            }

            return n_corr;
        }

    }  // namespace detail


    /**
     * Compute bedrock channel erosion during a single time step using the
     * Stream Power Law.
     *
     * It numerically solves the Stream Power Law [dh/dt = K A^m (dh/dx)^n]
     * using an implicit finite difference scheme 1st order in space and time.
     * The method is detailed in Braun and Willet's (2013) and has been
     * slightly reformulated/optimized.
     *
     * As it requires some input information on the flow tree topology
     * (e.g., flow receivers and stack order) and geometry (e.g., distance
     * to receivers), this generic function is grid/mesh agnostic and can
     * be applied in both 1-d (river profile) and 2-d (river network)
     * cases.
     *
     * The implicit scheme ensure numerical stability but doesn't totally
     * prevent erosion from lowering the elevation of a node below that of
     * its receiver. If this occurs, erosion will be limited so that the
     * node will be lowered (nearly) down to the level of its receiver. In
     * general it is better to adjust input values (e.g., ``dt``) so that
     * such arbitrary limitation doesn't occur. The value returned by this
     * function allows to detect and track the number of these
     * occurrences.
     *
     * @param erosion : ``[intent=out, shape=(nrows, ncols)||(nnodes)]``
     *     Erosion at grid node.
     * @param elevation : ``[intent=in, shape=(nrows, ncols)||(nnodes)]``
     *     Elevation at grid node.
     * @param flow_graph :``[intent=in]``
     *     Flow graph of the grid.
     * @param k_coef : ``[intent=in]``
     *     Stream Power Law coefficient.
     * @param m_exp : ``[intent=in]``
     *     Stream Power Law drainage area exponent.
     * @param n_exp : ``[intent=in]``
     *     Stream Power Law slope exponent.
     * @param dt : ``[intent=in]``
     *     Time step duration.
     * @param tolerance : ``[intent=in]``
     *     Tolerance used for Newton's iterations (``n_exp`` != 1).
     *
     * @returns
     *     Total number of nodes for which erosion has been
     *     arbitrarily limited to ensure consistency.
     */
    template <class Er, class El, class Dr, class FG>
    auto erode_stream_power(xtensor_t<Er>& erosion,
                            const xtensor_t<El>& elevation,
                            const xtensor_t<Dr>& drainage_area,
                            const FG& flow_graph,
                            double k_coef,
                            double m_exp,
                            double n_exp,
                            double dt,
                            double tolerance) -> typename FG::size_type
    {
        return detail::erode_stream_power_impl(erosion.derived_cast(),
                                               elevation.derived_cast(),
                                               drainage_area.derived_cast(),
                                               flow_graph,
                                               k_coef,
                                               m_exp,
                                               n_exp,
                                               dt,
                                               tolerance);
    }


    /**
     * Compute bedrock channel erosion during a single time step using the
     * Stream Power Law.
     *
     * This version accepts a spatially variable stream-power law coefficient.
     *
     * @param erosion : ``[intent=out, shape=(nrows, ncols)||(nnodes)]``
     *     Erosion at grid node.
     * @param elevation : ``[intent=in, shape=(nrows, ncols)||(nnodes)]``
     *     Elevation at grid node.
     * @param flow_graph :``[intent=in]``
     *     Flow graph of the grid.
     * @param k_coef : ``[intent=in, shape=(nrows, ncols)||(nnodes)]``
     *     Stream Power Law coefficient.
     * @param m_exp : ``[intent=in]``
     *     Stream Power Law drainage area exponent.
     * @param n_exp : ``[intent=in]``
     *     Stream Power Law slope exponent.
     * @param dt : ``[intent=in]``
     *     Time step duration.
     * @param tolerance : ``[intent=in]``
     *     Tolerance used for Newton's iterations (``n_exp`` != 1).
     *
     * @returns
     *     Total number of nodes for which erosion has been
     *     arbitrarily limited to ensure consistency.
     */
    template <class Er, class El, class Dr, class FG, class K>
    auto erode_stream_power(xtensor_t<Er>& erosion,
                            const xtensor_t<El>& elevation,
                            const xtensor_t<Dr>& drainage_area,
                            const FG& flow_graph,
                            const xtensor_t<K>& k_coef,
                            double m_exp,
                            double n_exp,
                            double dt,
                            double tolerance) -> typename FG::size_type
    {
        return detail::erode_stream_power_impl(erosion.derived_cast(),
                                               elevation.derived_cast(),
                                               drainage_area.derived_cast(),
                                               flow_graph,
                                               k_coef.derived_cast(),
                                               m_exp,
                                               n_exp,
                                               dt,
                                               tolerance);
    }

}  // namespace fastscapelib

#endif
