/**
 * xtensor utils (container tags, etc.)
 */
#ifndef FASTSCAPELIB_UTILS_XTENSOR_UTILS_H
#define FASTSCAPELIB_UTILS_XTENSOR_UTILS_H

#include <cstddef>

#include "xtensor/xarray.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xview.hpp"

namespace fastscapelib
{

    /**
     * The xtensor selector used by default in Fastscapelib C++ API.
     *
     * \rst
     *
     * Fastscapelib template classes specialized with this selector will use
     * :cpp:type:`xt::xtensor` or :cpp:type:`xt::xarray` as container types for
     * their public array members.
     *
     * \endrst
     *
     */
    struct xt_selector
    {
    };

    /**
     * Used to get the actual xtensor container type from a given selector.
     *
     * @tparam S The xtensor selector type.
     * @tparam T The container value type.
     * @tparam N The number of dimensions (only for static dimension containers)
     */
    template <class S, class T, std::size_t N = 0>
    struct container_selection
    {
    };

    template <class T, std::size_t N>
    struct container_selection<xt_selector, T, N>
    {
        using fixed_shape_type = xt::xtensor<T, N>;
        using dynamic_shape_type = xt::xarray<T>;
    };

    /**
     * Alias for the selected (static dimension) xtensor container type.
     *
     * @tparam S The container selector type.
     * @tparam T The container value type.
     * @tparam N The fixed number of dimensions.
     */
    template <class S, class T, std::size_t N>
    using fixed_shape_container_t = typename container_selection<S, T, N>::fixed_shape_type;

    /**
     * Alias for the selected (dynamic dimension) xtensor container type.
     *
     * @tparam S The xtensor selector type.
     * @tparam T The container value type.
     */
    template <class S, class T>
    using dynamic_shape_container_t = typename container_selection<S, T>::dynamic_shape_type;

    template <class T, class I>
    auto get_xt_view(T&& data, I&& max_idx)
    {
        return xt::view(xt::adapt(std::forward<T>(data)), xt::range(0, std::forward<I>(max_idx)));
    }
}  // namespace fastscapelib

#endif
