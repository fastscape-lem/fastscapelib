/**
 * Meta-programming stuff
 */
#pragma once

#include <cstdlib>

#include "xtensor/xtensor.hpp"
#include "xtensor/xarray.hpp"


namespace fastscapelib
{

/***********************
 * grid type selectors *
 ***********************/

struct profileG
{
};


struct rasterG
{
};


struct meshG
{
};


namespace detail
{


/*
 * Number of dimensions of any array that contains values at grid nodes.
 *
 * @tparam G Grid type selector.
 */
template <class G>
struct grid_nodes_ndim
{
    constexpr static std::size_t value = 1;
};

template <>
struct grid_nodes_ndim<rasterG>
{
    constexpr static std::size_t value = 2;
};


} // namespace detail


/*******************************
 * xtensor container selectors *
 *******************************/

struct xtensorC
{
};


struct xarrayC
{
};


namespace detail
{

template <class C, class T, std::size_t N = 0>
struct xt_container
{
};


template <class T, std::size_t N>
struct xt_container<xtensorC, T, N>
{
    using type = xt::xtensor<T, N>;
};


template <class T>
struct xt_container<xarrayC, T>
{
    using type = xt::xarray<T>;
};


} // namespace detail


template <class C, class T, std::size_t N = 0>
using xt_container_t = typename detail::xt_container<C, T, N>::type;


} // namespace fastscapelib
