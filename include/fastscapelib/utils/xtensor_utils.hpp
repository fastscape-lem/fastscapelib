/**
 * xtensor utils (container tags, etc.)
 */
#ifndef FASTSCAPELIB_UTILS_XTENSOR_UTILS_H
#define FASTSCAPELIB_UTILS_XTENSOR_UTILS_H


#include <cstddef>

#include "xtensor/xarray.hpp"
#include "xtensor/xtensor.hpp"


namespace fastscapelib
{

    struct xt_selector;

    template <class S, class T, std::size_t N = 0>
    struct xt_container;

    template <class T, std::size_t N>
    struct xt_container<xt_selector, T, N>
    {
        using tensor_type = xt::xtensor<T, N>;
        using array_type = xt::xarray<T>;
    };

    template <class S, class T, std::size_t N>
    using xt_tensor_t = typename xt_container<S, T, N>::tensor_type;

    template <class S, class T>
    using xt_array_t = typename xt_container<S, T>::array_type;

}  // namespace fastscapelib

#endif
