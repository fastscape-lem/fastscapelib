/**
 * Wrappers around RichDEM's routines.
 */
#pragma once

#ifdef ENABLE_RICHDEM

#include "richdem/richdem.hpp"


namespace fastscapelib
{

namespace detail
{

template<class A, class T = typename std::decay_t<A>::value_type>
auto to_array2d(A&& xt_array) -> richdem::Array2D<T>
{
    auto data = xt_array.data();
    auto shape = xt_array.shape();

    return richdem::Array2D<T>(data, shape[1], shape[0]);
}

}  // namespace detail


template<class E>
void fill_sinks_wei2018(xtensor_t<E>& elevation)
{
    auto elevation_a2d = detail::to_array2d(elevation.derived_cast());

    richdem::PriorityFlood_Wei2018(elevation_a2d);
}


}  // namespace fastscapelib

#endif
