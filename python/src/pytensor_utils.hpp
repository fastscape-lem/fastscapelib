#ifndef PYFASTSCAPELIB_PYTENSOR_UTILS_H
#define PYFASTSCAPELIB_PYTENSOR_UTILS_H

#include <cstddef>

#include "xtensor-python/pyarray.hpp"
#include "xtensor-python/pytensor.hpp"

#include "fastscapelib/utils/xtensor_utils.hpp"


namespace fastscapelib
{
    struct py_selector;

    template <class T, std::size_t N>
    struct xt_container<py_selector, T, N>
    {
        using tensor_type = xt::pytensor<T, N, xt::layout_type::row_major>;
        using array_type = xt::pyarray<T, xt::layout_type::row_major>;
    };
}

#endif
