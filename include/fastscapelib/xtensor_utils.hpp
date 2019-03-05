/**
 * xtensor utils (container tags, etc.)
 */
#pragma once

#include <cstddef>

#include "xtensor/xarray.hpp"
#include "xtensor/xtensor.hpp"


namespace fastscapelib
{

struct xtensor_tag;
struct xarray_tag;


template <class Tag, class T, std::size_t N = 0>
struct xt_container;


template <class T, std::size_t N>
struct xt_container<xtensor_tag, T, N>
{
    using type = xt::xtensor<T, N>;
};


template <class T>
struct xt_container<xarray_tag, T>
{
    using type = xt::xarray<T>;
};


template <class Tag, class T, std::size_t N = 0>
using xt_container_t = typename xt_container<Tag, T, N>::type;

} // namespace
