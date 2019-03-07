/**
 * Wrappers around RichDEM's routines.
 */
#pragma once

#include "richdem/common/Array2D.hpp"
#include "richdem/depressions/Wei2018.hpp"
#include "richdem/depressions/Zhou2016.hpp"
#include "richdem/flats/flats.hpp"
//#include "richdem/richdem.hpp"


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

template<class E>
void fill_sinks_zhou2016(xtensor_t<E>& elevation)
{
	auto elevation_a2d = detail::to_array2d(elevation.derived_cast());

	richdem::PriorityFlood_Zhou2016(elevation_a2d);
}


template<class E>
void resolve_flats_sloped(xtensor_t<E>& elevation)
{
    auto elevation_a2d = detail::to_array2d(elevation.derived_cast());

    richdem::ResolveFlatsEpsilon(elevation_a2d);
}


}  // namespace fastscapelib
