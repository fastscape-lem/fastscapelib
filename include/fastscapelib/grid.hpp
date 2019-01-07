#pragma once

#include <cstdlib>


namespace fastscapelib
{


/**
 * A selector for grid types.
 */
enum class grid_type
{
    profile = 0,
    raster,
    mesh
};


/*
 * Basic properties of a grid (depending on its type).
 */
template <grid_type G>
struct grid_properties
{
    static constexpr std::size_t ndim_nodes = 1;
};


template <>
struct grid_properties<grid_type::raster>
{
    static constexpr std::size_t ndim_nodes = 2;
};


}  // namespace fastscapelib
