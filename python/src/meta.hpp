/**
 * Meta-programming stuff
 */
#pragma once

#include <cstdlib>

#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyarray.hpp"


namespace fastscapelib
{


/**************************************
 * xtensor-python container selectors *
 **************************************/

struct pytensorC
{
};


struct pyarrayC
{
};


namespace detail
{


template <class T, std::size_t N>
struct xt_container<pytensorC, T, N>
{
    using type = xt::pytensor<T, N>;
};


template <class T>
struct xt_container<pyarrayC, T>
{
    using type = xt::pyarray<T>;
};


} // namespace detail

} // namespace fastscapelib
