/**
 * @file
 * @brief Some utilities used internally.
 */
#pragma once

#include <cstddef>
#include <type_traits>

#include "xtensor/xcontainer.hpp"


// type used for array indexers
using index_t = std::ptrdiff_t;

// type used when working with both shapes or sizes (unsigned) and indexers (signed)
using sshape_t = std::make_signed_t<size_t>;

// type used as an alias to xtensor's xexpression
// TODO: use xexpression shaped when it's ready (see xtensor GH #994)
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
