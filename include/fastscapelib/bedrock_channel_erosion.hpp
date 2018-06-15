/**
 * Functions to compute bedrock channel erosion.
 */
#pragma once

#include <cmath>
#include <limits>

#include "xtensor/xtensor.hpp"
#include "xtensor/xstrided_view.hpp"


namespace fastscapelib
{

namespace detail
{


/**
 * erode_stream_power implementation.
 */
template<class Er, class El, class S, class R, class Di, class Dr>
void erode_stream_power_impl(Er&& erosion,
                             El&& elevation,
                             S&& stack,
                             R&& receivers,
                             Di&& dist2receivers,
                             Dr&& drainage_area,
                             double k_coef,
                             double m_exp,
                             double n_exp,
                             double dt,
                             double tolerance)
{
    using T = typename std::decay_t<El>::value_type;;

    auto erosion_flat = xt::flatten(erosion);
    const auto elevation_flat = xt::flatten(elevation);
    const auto drainage_area_flat = xt::flatten(drainage_area);

    for (auto&& istack : stack)
    {
        const index_t irec = receivers(istack);

        if (irec == istack)
        {
            // no erosion at basin outlets
            erosion(istack) = 0.0;
            continue;
        }

        const double factor = (k_coef * dt *
                               std::pow(drainage_area_flat(istack), m_exp) /
                               std::pow(dist2receivers(istack), n_exp));

        const T istack_elevation = elevation_flat(istack);
        const T irec_elevation = elevation_flat(irec) - erosion(irec);

        // iterate: lower elevation until convergence
        T elevation_k = istack_elevation;
        T elevation_prev = std::numeric_limits<T>::max();

        while (std::abs(elevation_k - elevation_prev) > tolerance)
        {
            const auto slope = elevation_k - irec_elevation;
            const auto diff = (
                (elevation_k - istack_elevation + factor * std::pow(slope, n_exp)) /
                (1. + factor * n_exp * std::pow(slope, n_exp - 1)));

            elevation_k -= diff;
            elevation_prev = elevation_k;
        }

        erosion_flat(istack) = istack_elevation - elevation_k;
    }
}

}  // namespace detail


/**
 * Compute bedrock channel erosion during a single time step using the
 * Stream Power Law.
 *
 * It uses an implicit scheme detailed in Braun and Willet's (2013).
 *
 * It supports both unstructured mesh and 2-d grid.
 *
 * @param erosion : ``[intent=out, shape=(nrows, ncols)||(nnodes)]``
 *     Erosion at grid node.
 * @param elevation : ``[intent=out, shape=(nrows, ncols)||(nnodes)]``
 *     Elevation at grid node.
 * @param stack :``[intent=in, shape=(nnodes)]``
 *     Stack position at grid node.
 * @param receivers : ``[intent=in, shape=(nnodes)]``
 *     Index of flow receiver at grid node.
 * @param dist2receivers : ``[intent=out, shape=(nnodes)]``
 *     Distance to receiver at grid node.
 * @param drainage_area : ``[intent=out, shape=(nrows, ncols)||(nnodes)]``
 *     Drainage area at grid node.
 * @param k_coef : ``[intent=in]``
 *     Stream Power Law coefficient.
 * @param m_exp : ``[intent=in]``
 *     Stream Power Law drainage area exponent.
 * @param n_exp : ``[intent=in]``
 *     Stream Power Law drainage slope exponent.
 * @param dt : ``[intent=in]``
 *     Time step duration.
 * @param tolerance : ``[intent=in]``
 *     Tolerance used for convergence of Newton's iterations.
 */
template<class Er, class El, class S, class R, class Di, class Dr>
void erode_stream_power(xtensor_t<Er>& erosion,
                        xtensor_t<El>& elevation,
                        xtensor_t<S>& stack,
                        xtensor_t<R>& receivers,
                        xtensor_t<Di>& dist2receivers,
                        xtensor_t<Dr>& drainage_area,
                        double k_coef,
                        double m_exp,
                        double n_exp,
                        double dt,
                        double tolerance)
{
    detail::erode_stream_power_impl(erosion.derived_cast(),
                                    elevation.derived_cast(),
                                    stack.derived_cast(),
                                    receivers.derived_cast(),
                                    dist2receivers.derived_cast(),
                                    drainage_area.derived_cast(),
                                    k_coef, m_exp, n_exp,
                                    dt, tolerance);
}

}  // namespace fastscapelib
