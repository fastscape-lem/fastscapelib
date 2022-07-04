/**
 * Functions to compute bedrock channel erosion.
 */
#ifndef FASTSCAPELIB_ERODERS_SPL_H
#define FASTSCAPELIB_ERODERS_SPL_H

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <type_traits>

#include "xtensor/xbroadcast.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xmanipulation.hpp"

#include "fastscapelib/utils/utils.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"


namespace fastscapelib
{

    /**
     * Bedrock channel erosion modelled using the Stream Power Law.
     *
     * It numerically solves the Stream Power Law [dh/dt = K A^m (dh/dx)^n]
     * using an implicit finite difference scheme 1st order in space and time.
     * The method is detailed in Braun and Willet's (2013) and has been
     * slightly reformulated/optimized.
     *
     * The implicit scheme ensure numerical stability but doesn't totally
     * prevent erosion from lowering the elevation of a node below that of its
     * receiver. If this occurs, erosion will be limited so that the node will
     * be lowered (nearly) down to the level of its receiver. In general it is
     * better to adjust input values (e.g., time step or input parameters) so
     * that such arbitrary limitation doesn't occur.
     *
     * `n_corr` is the total number of nodes for which erosion has been
     * arbitrarily limited to ensure consistency.
     */
    template <class FG, class S = typename FG::xt_selector>
    class spl_eroder
    {
    public:
        using flow_graph_type = FG;
        using xt_selector = S;

        using size_type = typename flow_graph_type::size_type;
        using shape_type = typename flow_graph_type::shape_type;
        using data_type = typename flow_graph_type::data_type;
        using data_array_type = xt_array_t<xt_selector, data_type>;

        template <class K>
        spl_eroder(FG& flow_graph, K&& k_coef, double area_exp, double slope_exp, double tolerance)
            : m_flow_graph(flow_graph)
            , m_shape(flow_graph.grid_shape())
            , m_tolerance(tolerance)
        {
            set_k_coef(k_coef);
            set_area_exp(area_exp);
            set_slope_exp(slope_exp);

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
            if (!xt::same_shape(value.shape(), m_shape))
            {
                throw std::runtime_error("cannot set k_coef value: shape mismatch");
            }
            m_k_coef = xt::flatten(value);
        };

        double area_exp()
        {
            return m_area_exp;
        };

        void set_area_exp(double value)
        {
            // TODO: validate value
            m_area_exp = value;
        };

        double slope_exp()
        {
            return m_slope_exp;
        };

        void set_slope_exp(double value)
        {
            // TODO: validate value
            m_slope_exp = value;
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
        double m_area_exp;
        double m_slope_exp;
        double m_tolerance;
        size_type m_n_corr;
    };

    /**
     * SPL erosion implementation.
     *
     * The implementation here slightly differs from the one described in
     * Braun & Willet (2013). The problem is reformulated so that the
     * Newton-Raphson method is applied on the difference of elevation
     * between a node and its receiver, rather than on the node's
     * elevation itself. This allows saving some operations.
     */
    template <class FG, class S>
    auto spl_eroder<FG, S>::erode(const data_array_type& elevation,
                                  const data_array_type& drainage_area,
                                  double dt) -> const data_array_type&
    {
        auto& flow_graph_impl = m_flow_graph.impl();

        const auto& receivers = flow_graph_impl.receivers();
        const auto& receivers_count = flow_graph_impl.receivers_count();
        const auto& receivers_distance = flow_graph_impl.receivers_distance();
        const auto& receivers_weight = flow_graph_impl.receivers_weight();
        const auto& dfs_indices = flow_graph_impl.dfs_indices();

        m_n_corr = 0;

        for (const auto& idfs : dfs_indices)
        {
            auto r_count = receivers_count[idfs];

            for (size_type r = 0; r < r_count; ++r)
            {
                size_type irec = receivers(idfs, r);

                if (irec == idfs)
                {
                    // at basin outlet or pit
                    m_erosion.flat(idfs) = 0.;
                    continue;
                }

                data_type idfs_elevation = elevation.flat(idfs);  // at time t
                data_type irec_elevation
                    = elevation.flat(irec) - m_erosion.flat(irec);  // at time t+dt

                if (irec_elevation >= idfs_elevation)
                {
                    // may happen if flow is routed outside of a depression / flat area
                    m_erosion.flat(idfs) = 0.;
                    continue;
                }

                auto factor
                    = (m_k_coef(idfs) * dt * std::pow(drainage_area.flat(idfs), m_area_exp));

                data_type delta_0 = idfs_elevation - irec_elevation;
                data_type delta_k;

                if (m_slope_exp == 1)
                {
                    // fast path for slope_exp = 1 (common use case)
                    factor /= receivers_distance(idfs, r);
                    delta_k = delta_0 / (1. + factor);
                }

                else
                {
                    // 1st order Newton-Raphson iterations (k)
                    factor /= std::pow(receivers_distance(idfs, r), m_slope_exp);
                    delta_k = delta_0;

                    // TODO: add convergence control parameters (max_iterations, atol, mtol)
                    while (true)
                    {
                        auto factor_delta_exp = factor * std::pow(delta_k, m_slope_exp);
                        auto func = delta_k + factor_delta_exp - delta_0;

                        if (func <= m_tolerance)
                        {
                            break;
                        }

                        auto func_deriv = 1 + m_slope_exp * factor_delta_exp / delta_k;
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
                    m_erosion.flat(idfs) = (delta_0 - delta_k);
                }
                else
                {
                    m_erosion.flat(idfs) += (delta_0 - delta_k) * receivers_weight(idfs, r);
                }
            }
        }

        return m_erosion;
    };


    template <class FG, class K>
    spl_eroder<FG> make_spl_eroder(
        FG& graph, K&& k_coef, double area_exp, double slope_exp, double tolerance)
    {
        return spl_eroder<FG>(graph, k_coef, area_exp, slope_exp, tolerance);
    }
}  // namespace fastscapelib

#endif
