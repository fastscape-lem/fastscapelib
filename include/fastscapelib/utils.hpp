/**
 * @file
 * @brief Some utilities used internally.
 */
#pragma once

#include <cstddef>

#include "xtensor/xcontainer.hpp"


// type used for indexing arrays and array sizes.
// TODO: use xt::index_t directly if/when xtensor #661 is merged.
using index_t = std::ptrdiff_t;


template<class E>
using xtensor_t = xt::xexpression<E>;


namespace fastscapelib
{

namespace detail
{


/**
 * @brief Return true if a given (row, col) index is in bounds (false
 *        otherwise).
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


}  // namespace detail

}  // namespace fastscapelib
