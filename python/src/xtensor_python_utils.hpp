#pragma once

#include <cstdlib>

#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyarray.hpp"

#include "fastscapelib/xtensor_utils.hpp"


namespace fastscapelib
{


template <class T, std::size_t N>
struct xt_container<array_type::pytensor, T, N>
{
    using type = xt::pytensor<T, N>;
};


template <class T>
struct xt_container<array_type::pyarray, T>
{
    using type = xt::pyarray<T>;
};


} // namespace fastscapelib
