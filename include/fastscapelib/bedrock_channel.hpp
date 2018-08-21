/**
 * Functions to compute bedrock channel erosion.
 */
#pragma once

#include <cassert>
#include <cmath>
#include <limits>
#include <type_traits>

#include "xtensor/xtensor.hpp"
#include "xtensor/xstrided_view.hpp"

#include "fastscapelib/utils.hpp"


namespace fastscapelib
{

namespace detail
{


template<class T, class S>
auto k_coef_as_array(T k_coef,
                     S&& shape,
                     typename std::enable_if_t<std::is_floating_point<T>::value>* = 0)
{
    return xt::broadcast(std::forward<T>(k_coef), shape);
}


template<class K, class S>
auto k_coef_as_array(K&& k_coef,
                     S&& shape,
                     typename std::enable_if_t<xt::is_xexpression<K>::value>* = 0)
{
    auto k_coef_arr = xt::flatten(k_coef);

    assert(k_coef_arr.shape() == shape);
    (void) shape;   // TODO: still unused parameter warning despite assert?

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
template<class Er, class El, class S, class R, class Di, class Dr, class K>
index_t erode_stream_power_impl(Er&& erosion,
                                El&& elevation,
                                S&& stack,
                                R&& receivers,
                                Di&& dist2receivers,
                                Dr&& drainage_area,
                                K&& k_coef,
                                double m_exp,
                                double n_exp,
                                double dt,
                                double tolerance)
{
    using T = std::common_type_t<typename std::decay_t<Er>::value_type,
                                 typename std::decay_t<El>::value_type>;

    auto erosion_flat = xt::flatten(erosion);
    const auto elevation_flat = xt::flatten(elevation);
    const auto drainage_area_flat = xt::flatten(drainage_area);

    const auto k_coef_arr = k_coef_as_array(k_coef, elevation_flat.shape());

    index_t n_corr = 0;

    for (const auto& istack : stack)
    {
        index_t irec = receivers(istack);

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

        auto factor = (k_coef_arr(istack) * dt *
                       std::pow(drainage_area_flat(istack), m_exp));

        T delta_0 = istack_elevation - irec_elevation;
        T delta_k;

        if (n_exp == 1)
        {
            // fast path for n_exp = 1 (common use case)
            factor /= dist2receivers(istack);
            delta_k = delta_0 / (1. + factor);
        }

        else
        {
            // 1st order Newton-Raphson iterations (k)
            factor /= std::pow(dist2receivers(istack), n_exp);
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

        erosion_flat(istack) = delta_0 - delta_k;
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
template<class Er, class El, class S, class R, class Di, class Dr>
index_t erode_stream_power(xtensor_t<Er>& erosion,
                           const xtensor_t<El>& elevation,
                           const xtensor_t<S>& stack,
                           const xtensor_t<R>& receivers,
                           const xtensor_t<Di>& dist2receivers,
                           const xtensor_t<Dr>& drainage_area,
                           double k_coef,
                           double m_exp,
                           double n_exp,
                           double dt,
                           double tolerance)
{
    return detail::erode_stream_power_impl(erosion.derived_cast(),
                                           elevation.derived_cast(),
                                           stack.derived_cast(),
                                           receivers.derived_cast(),
                                           dist2receivers.derived_cast(),
                                           drainage_area.derived_cast(),
                                           k_coef, m_exp, n_exp,
                                           dt, tolerance);
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
 * @param stack :``[intent=in, shape=(nnodes)]``
 *     Stack position at grid node.
 * @param receivers : ``[intent=in, shape=(nnodes)]``
 *     Index of flow receiver at grid node.
 * @param dist2receivers : ``[intent=out, shape=(nnodes)]``
 *     Distance to receiver at grid node.
 * @param drainage_area : ``[intent=out, shape=(nrows, ncols)||(nnodes)]``
 *     Drainage area at grid node.
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
template<class Er, class El, class S, class R, class Di, class Dr, class K>
index_t erode_stream_power(xtensor_t<Er>& erosion,
                           const xtensor_t<El>& elevation,
                           const xtensor_t<S>& stack,
                           const xtensor_t<R>& receivers,
                           const xtensor_t<Di>& dist2receivers,
                           const xtensor_t<Dr>& drainage_area,
                           const xtensor_t<K>& k_coef,
                           double m_exp,
                           double n_exp,
                           double dt,
                           double tolerance)
{
    return detail::erode_stream_power_impl(erosion.derived_cast(),
                                           elevation.derived_cast(),
                                           stack.derived_cast(),
                                           receivers.derived_cast(),
                                           dist2receivers.derived_cast(),
                                           drainage_area.derived_cast(),
                                           k_coef.derived_cast(),
                                           m_exp, n_exp,
                                           dt, tolerance);
}

}  // namespace fastscapelib
