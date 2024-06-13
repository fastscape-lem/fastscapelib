#ifndef PYFASTSCAPELIB_PYTENSOR_CONTAINERS_HPP
#define PYFASTSCAPELIB_PYTENSOR_CONTAINERS_HPP

#include "fastscapelib/utils/containers.hpp"

#include "xtensor-python/pyarray.hpp"
#include "xtensor-python/pytensor.hpp"

#include <cstddef>


namespace fastscapelib
{
    struct xt_python_selector;

    template <class T, std::size_t N>
    struct container_selection<xt_python_selector, T, N>
    {
        using fixed_shape_type = xt::pytensor<T, N, xt::layout_type::row_major>;
        using dynamic_shape_type = xt::pyarray<T, xt::layout_type::row_major>;
    };

    template <class T, xt::layout_type L>
    struct container_impl<xt::pyarray<T, L>> : public xtensor_shared_utils<xt::pyarray<T, L>>
    {
        using base_type = xtensor_shared_utils<xt::pyarray<T, L>>;
        using container_type = typename xt::pyarray<T, L>;
        using size_type = typename container_type::size_type;
        using shape_type = typename container_type::shape_type;

        using base_type::compute_distance;
        using base_type::get_view;
        using base_type::init;
    };

    template <class T, std::size_t N, xt::layout_type L>
    struct container_impl<xt::pytensor<T, N, L>>
        : public xtensor_shared_utils<xt::pytensor<T, N, L>>
    {
        using base_type = xtensor_shared_utils<xt::pytensor<T, N, L>>;
        using container_type = typename xt::pytensor<T, N, L>;
        using size_type = typename container_type::size_type;
        using shape_type = typename container_type::shape_type;

        using base_type::compute_distance;
        using base_type::get_view;
        using base_type::init;
    };
}

#endif  // PYFASTSCAPELIB_PYTENSOR_CONTAINERS_HPP
