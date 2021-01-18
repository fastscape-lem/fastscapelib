/**
 * xtensor utils (container tags, etc.)
 */
#ifndef FASTSCAPELIB_XTENSOR_UTILS_H
#define FASTSCAPELIB_XTENSOR_UTILS_H


#include <cstddef>

#include "xtensor/xarray.hpp"
#include "xtensor/xtensor.hpp"


namespace fastscapelib
{

    struct xtensor_selector;
    struct xarray_selector;


    template <class X, class T, std::size_t N = 0>
    struct xt_container;


    template <class T, std::size_t N>
    struct xt_container<xtensor_selector, T, N>
    {
        using type = xt::xtensor<T, N>;
    };


    template <class T>
    struct xt_container<xarray_selector, T>
    {
        using type = xt::xarray<T>;
    };


    template <class X, class T, std::size_t N = 0>
    using xt_container_t = typename xt_container<X, T, N>::type;

} // namespace fastscapelib

#endif