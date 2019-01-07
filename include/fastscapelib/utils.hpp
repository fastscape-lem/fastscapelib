/**
 * @file
 * @brief Some utilities used internally.
 */
#pragma once

#include <cstddef>
#include <type_traits>

#include "xtensor/xcontainer.hpp"

#include "fastscapelib/xtensor_utils.hpp"


namespace fastscapelib
{


namespace utils
{


/**
 * Return true if a given (row, col) index is in bounds (false otherwise).
 */
template<class S>
bool in_bounds(const S& shape, index_t row, index_t col)
{
    if(row >= 0 && row < static_cast<index_t>(shape[0])
       && col >= 0 && col < static_cast<index_t>(shape[1]))
    {
        return true;
    }
    return false;
}


/**
 * modulus (note: k % n is division remainder).
 */
template<class T>
T mod(T k, T n)
{
    return ((k %= n) < 0) ? k+n : k;
}


/**
 * Static cast to strongly typed enum's underlying type.
 */
template <typename E>
constexpr auto cast_underlying(E e) noexcept
{
    return static_cast<std::underlying_type_t<E>>(e);
}


}  // namespace utils


}  // namespace fastscapelib
