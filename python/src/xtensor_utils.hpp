/**
 * xtensor utils (container tags, etc.)
 */
#pragma once

#include <cstddef>

#include "xtensor-python/pyarray.hpp"
#include "xtensor-python/pytensor.hpp"

#include "fastscapelib/xtensor_utils.hpp"


namespace fastscapelib
{

struct pytensor_selector;
struct pyarray_selector;

template <class T, std::size_t N>
struct xt_container<pytensor_selector, T, N>
{
    using type = xt::pytensor<T, N>;
};

template <class T>
struct xt_container<pyarray_selector, T>
{
    using type = xt::pyarray<T>;
};

} // namespace fastscapelib
