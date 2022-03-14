#ifndef PYFASTSCAPELIB_PYTENSOR_UTILS_H
#define PYFASTSCAPELIB_PYTENSOR_UTILS_H

#include <cstddef>

#include "xtensor-python/pyarray.hpp"

#include "fastscapelib/xtensor_utils.hpp"


namespace fastscapelib
{
    struct pyarray_selector;

    template <class T, std::size_t N>
    struct xt_container<pyarray_selector, T, N>
    {
        using type = xt::pyarray<T, xt::layout_type::row_major>;
    };
}

#endif
