
#ifndef FASTSCAPELIB_TEST_RASTER_GRID_H
#define FASTSCAPELIB_TEST_RASTER_GRID_H

#include "fastscapelib/grid/raster_grid.hpp"

#include "gtest/gtest.h"

#include <array>


namespace fs = fastscapelib;


namespace fastscapelib
{
    namespace testing
    {

        class raster_boundary_status : public ::testing::Test
        {
        protected:
            using node_s = fs::node_status;

            node_status fb = node_status::fixed_value;
            node_status lb = node_status::looped;
            node_status co = node_status::core;

            std::array<node_s, 4> hloop{ { lb, lb, fb, fb } };
            std::array<node_s, 4> vloop{ { fb, fb, lb, lb } };
            std::array<node_s, 4> hvloop{ { lb, lb, lb, lb } };
            std::array<node_s, 4> hill_formed_loop{ { lb, fb, lb, lb } };

            fs::raster_boundary_status fixed_value_status{ fb };
            fs::raster_boundary_status hlooped_status{ hloop };
            fs::raster_boundary_status vlooped_status{ vloop };
            fs::raster_boundary_status hvlooped_status{ hvloop };
        };

        class raster_grid_base : public raster_boundary_status
        {
        protected:
            using size_type = typename fs::raster_grid::size_type;
            using length_type = typename fs::raster_grid::length_type;
            using spacing_type = typename fs::raster_grid::spacing_type;
            using shape_type = typename fs::raster_grid::shape_type;

            shape_type shape{ { 5, 10 } };
        };

        class raster_grid : public raster_grid_base
        {
        protected:
            using grid_type = fs::raster_grid;

            grid_type fixed_grid = grid_type(shape, { 1.3, 1.2 }, fixed_value_status);
            grid_type hlooped_grid = grid_type(shape, { 1.3, 1.2 }, hlooped_status);
            grid_type vlooped_grid = grid_type(shape, { 1.3, 1.2 }, vlooped_status);
            grid_type looped_grid = grid_type(shape, { 1.4, 1.8 }, hvlooped_status);
        };

        class queen_raster_grid : public raster_grid_base
        {
        protected:
            using grid_type = fs::raster_grid_xt<fs::xt_selector, fs::raster_connect::queen>;

            grid_type queen_fixed = grid_type(shape, { 1.3, 1.2 }, fixed_value_status);
            grid_type queen_hlooped = grid_type(shape, { 1.3, 1.2 }, hlooped_status);
            grid_type queen_vlooped = grid_type(shape, { 1.3, 1.2 }, vlooped_status);
            grid_type queen_looped = grid_type(shape, { 1.4, 1.8 }, hvlooped_status);

            grid_type::neighbors_indices_type neighbors_idx;
        };

        class queen_raster_grid_fixed : public queen_raster_grid
        {
        };

        class queen_raster_grid_looped : public queen_raster_grid
        {
        };

        class rook_raster_grid : public raster_grid_base
        {
        protected:
            using grid_type = fs::raster_grid_xt<fs::xt_selector, fs::raster_connect::rook>;

            grid_type rook_fixed = grid_type(shape, { 1.3, 1.2 }, fixed_value_status);
            grid_type rook_hlooped = grid_type(shape, { 1.3, 1.2 }, hlooped_status);
            grid_type rook_vlooped = grid_type(shape, { 1.3, 1.2 }, vlooped_status);
            grid_type rook_looped = grid_type(shape, { 1.4, 1.8 }, hvlooped_status);

            grid_type::neighbors_indices_type neighbors_idx;
        };

        class rook_raster_grid_fixed : public rook_raster_grid
        {
        };

        class rook_raster_grid_looped : public rook_raster_grid
        {
        };

        class bishop_raster_grid : public raster_grid_base
        {
        protected:
            using grid_type = fs::raster_grid_xt<fs::xt_selector, fs::raster_connect::bishop>;

            grid_type bishop_fixed = grid_type(shape, { 1.3, 1.2 }, fixed_value_status);
            grid_type bishop_hlooped = grid_type(shape, { 1.3, 1.2 }, hlooped_status);
            grid_type bishop_vlooped = grid_type(shape, { 1.3, 1.2 }, vlooped_status);
            grid_type bishop_looped = grid_type(shape, { 1.4, 1.8 }, hvlooped_status);

            grid_type::neighbors_indices_type neighbors_idx;
        };

        class bishop_raster_grid_fixed : public bishop_raster_grid
        {
        };

        class bishop_raster_grid_looped : public bishop_raster_grid
        {
        };
    }
}

#endif
