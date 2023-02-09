/**
 * Functions to compute bedrock channel erosion.
 */
#ifndef FASTSCAPELIB_ERODERS_SPL_H
#define FASTSCAPELIB_ERODERS_SPL_H

#include <algorithm>
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
     * slightly adapted.
     *
     * For the linear case (n = 1), solving the equation is trivial. For the
     * non-linear case (n != 1), the Newton-Raphson method is used to find the
     * optimal solution.
     *
     * This implementation supports multiple direction flow, although only for
     * the linear case.
     *
     * It also prevents the formation of new closed depressions (lakes). The
     * erosion within existing closed depressions of the input topographic
     * surface will depend on the given flow graph:
     *
     * - If the flow graph handles routing accross those depressions, no erosion
     *   will take place within them (lakes).
     * - If the flow graph does not resolve those depressions (no sink resolver
     *   applied), they will be eroded like the rest of the topographic surface
     *   (no lake).
     *
     * Caveat: the numerical scheme used here is stable but not implicit with
     * respect to upslope contributing area (A). Using large time steps may
     * still have a great impact on the solution, especially on transient
     * states.
     *
     * @tparam FG The flow graph type.
     * @tparam S The xtensor container selector for data array members.
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
        spl_eroder(
            FG& flow_graph, K&& k_coef, double area_exp, double slope_exp, double tolerance = 1e-3)
            : m_flow_graph(flow_graph)
            , m_shape(flow_graph.grid_shape())
            , m_tolerance(tolerance)
        {
            set_k_coef(k_coef);
            set_area_exp(area_exp);
            set_slope_exp(slope_exp);

            m_erosion.resize(m_shape);
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
            // TODO: do not allow n != 1 with multiple flow
            m_slope_exp = value;
            m_linear = (std::fabs(value) - 1) <= std::numeric_limits<double>::epsilon();
        };

        double tolerance()
        {
            return m_tolerance;
        }

        /**
         * Returns the number of nodes for which erosion has been arbitrarily
         * limited during the last computation.
         *
         * To ensure numerical stability, channel erosion may not lower the
         * elevation of a node below so that it reverts the slope with its
         * direct neighbors. This prevents the formation of new closed
         * depressions.
         */
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
        bool m_linear;
        size_type m_n_corr;
    };

    /**
     * SPL erosion implementation.
     *
     * The implementation here slightly differs from the one described in Braun
     * & Willet (2013). For the non-linear case, the problem is reformulated so
     * that the Newton-Raphson method is applied on the difference of elevation
     * between a node and its receiver, rather than on the node's elevation
     * itself.
     *
     * This implementation also works when the input elevation is not fully
     * consistent with flow paths: if the topographic surface still contains
     * closed depressions that have been resolved in the current flow graph state,
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

        // reset
        m_erosion.fill(0);
        m_n_corr = 0;

        // iterate over graph nodes in the bottom->up direction
        for (const auto& inode : flow_graph_impl.node_indices_bottomup())
        {
            data_type inode_elevation = elevation.flat(inode);
            auto r_count = receivers_count[inode];

            if (r_count == 1 && receivers(inode, 0) == inode)
            {
                // current node is a basin outlet or pit (no erosion)
                continue;
            }

            // ``elevation_flooded`` represents the level below which no erosion
            // should happen (the node is already within a lake or in order
            // prevent the formation of new lakes). It corresponds to the lowest
            // elevation at step + 1 found among the node receivers
            double elevation_flooded = std::numeric_limits<double>::max();

            for (size_type r = 0; r < r_count; ++r)
            {
                size_type irec = receivers(inode, r);
                data_type irec_elevation_next = elevation.flat(irec) - m_erosion.flat(irec);

                if (irec_elevation_next < elevation_flooded)
                {
                    elevation_flooded = irec_elevation_next;
                }
            }

            if (inode_elevation <= elevation_flooded)
            {
                // current node is inside a lake (no erosion)
                continue;
            }

            // init discrete equation numerator and denominator values
            double eq_num = inode_elevation;
            double eq_den = 1.0;

            // iterate over receivers
            for (size_type r = 0; r < r_count; ++r)
            {
                size_type irec = receivers(inode, r);

                data_type irec_elevation = elevation.flat(irec);
                // elevation of receiver node at step + 1
                data_type irec_elevation_next = irec_elevation - m_erosion.flat(irec);

                if (irec_elevation > inode_elevation)
                {
                    // current node is next to a lake spill ;
                    // the current receiver is not the lake spill, it doesn't contribute to
                    // eroding the current node
                    continue;
                }

                data_type irec_weight = receivers_weight(inode, r);
                data_type irec_distance = receivers_distance(inode, r);

                auto factor = (m_k_coef(inode) * dt
                               * std::pow(drainage_area.flat(inode) * irec_weight, m_area_exp));

                if (m_linear)
                {
                    // fast path for the linear case
                    factor /= irec_distance;
                    eq_num += factor * irec_elevation_next;
                    eq_den += factor;
                }
                else
                {
                    // 1st order Newton-Raphson iterations for the non-linear case
                    factor /= std::pow(irec_distance, m_slope_exp);

                    // solve directly for the difference of elevation
                    // faster?
                    // (only valid for single flow direction)
                    double delta_0 = inode_elevation - irec_elevation_next;
                    double delta_k = delta_0;

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

                    eq_num = inode_elevation - (delta_0 - delta_k);
                }
            }

            data_type inode_elevation_updated = eq_num / eq_den;

            if (inode_elevation_updated < elevation_flooded)
            {
                // numerical stability: prevent the creation of new
                // depressions / flat channels by arbitrarily limiting
                // erosion
                m_n_corr++;
                inode_elevation_updated = elevation_flooded + std::numeric_limits<data_type>::min();
            }

            m_erosion.flat(inode) = inode_elevation - inode_elevation_updated;
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
