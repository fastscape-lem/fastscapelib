#ifndef PYFASTSCAPELIB_PYTENSOR_UTILS_H
#define PYFASTSCAPELIB_PYTENSOR_UTILS_H

#include <cstddef>


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
