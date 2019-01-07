#pragma once

#include <cstddef>
#include <type_traits>

#include "xtensor/xcontainer.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xarray.hpp"


namespace fastscapelib
{


/**
 * A selector for xtensor's array (container) types.
 */
enum class array_type
{
    xtensor = 0,
    xarray,
    pytensor,
    pyarray
};


template <array_type A, class T, std::size_t N = 0>
struct xt_container
{
};


template <class T, std::size_t N>
struct xt_container<array_type::xtensor, T, N>
{
    using type = xt::xtensor<T, N>;
};


template <class T>
struct xt_container<array_type::xarray, T>
{
    using type = xt::xarray<T>;
};


template <array_type A, class T, std::size_t N = 0>
using xt_container_t = typename xt_container<A, T, N>::type;


/**
 * Type used for array indexers.
 */
using index_t = std::ptrdiff_t;


/**
 * Type used when working with both shapes or sizes (unsigned) and indexers (signed).
 */
using sshape_t = std::make_signed_t<std::size_t>;


/**
 * Alias of xtensor's xexpression, used just to make function signatures a little
 * more explicit.
 *
 * TODO: use xexpression shaped when it's ready (see xtensor GH #994)
 */
template<class E>
using xtensor_t = xt::xexpression<E>;


}  // namespace fastscapelib
