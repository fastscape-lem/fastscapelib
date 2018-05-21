/**
 * @file
 * @brief Some utilities used internally.
 */
#pragma once

#include <cstdint>

#include "xtensor/xcontainer.hpp"


// type used for indexing arrays and array sizes.
using index_t = int64_t;


template<typename D>
using xtensor_t = xt::xcontainer<D>;


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
    if(row >= 0 && row < (index_t) shape[0]
       && col >= 0 && col < (index_t) shape[1])
    {
        return true;
    }
    return false;
}


}  // namespace detail

}  // namespace fastscapelib
