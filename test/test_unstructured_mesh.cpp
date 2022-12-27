#include <stdexcept>

#include "gtest/gtest.h"
#include "xtensor/xindex_view.hpp"
#include "xtensor/xtensor.hpp"

#include "fastscapelib/grid/unstructured_mesh.hpp"


namespace fs = fastscapelib;


namespace fastscapelib
{
    namespace testing
    {

        class unstructured_mesh : public ::testing::Test
        {
        protected:
            using node_s = fs::node_status;
            node_s fixed = node_s::fixed_value_boundary;

            using grid_type = fs::unstructured_mesh_xt<fs::xt_selector>;
            using size_type = typename grid_type::size_type;
            using shape_type = typename grid_type::shape_type;

            size_type mesh_size = 5;

            // See Python tests for details about the mesh arrays used here
            xt::xtensor<double, 2> points{
                { 0.0, 0.5 }, { 0.5, 0.0 }, { 0.0, -0.5 }, { -0.5, 0.0 }, { 0.0, 0.0 }
            };

            xt::xtensor<size_type, 1> indptr{ 0, 3, 6, 9, 12, 16 };
            xt::xtensor<size_type, 1> indices{ 4, 3, 1, 4, 2, 0, 4, 3, 1, 2, 4, 0, 2, 3, 1, 0 };
            xt::xtensor<size_type, 1> convex_hull_indices{ 0, 1, 2, 3 };

            xt::xtensor<double, 1> areas{ 1.0, 1.0, 1.0, 1.0, 2.0 };

            xt::xtensor<double, 2> neighbors_distances
            {
                { 0.5, 0.70710678, 0.70710678 0. }, { 0.5, 0.70710678, 0.70710678, 0. },
                    { 0.5, 0.70710678, 0.70710678, 0. }, { 0.70710678, 0.5, 0.70710678, 0. },
                {
                    0.5, 0.5, 0.5, 0.5
                }
            }
        }
    };

    TEST_F(unstructured_mesh, static_expr)
    {
        EXPECT_EQ(fs::unstructured_mesh::is_structured(), false);
        EXPECT_EQ(fs::unstructured_mesh::is_uniform(), false);
        EXPECT_EQ(fs::unstructured_mesh::n_neighbors_max(), 30u);
        EXPECT_EQ(fs::unstructured_mesh::xt_ndims(), 1);
    }


    TEST_F(unstructured_mesh, ctor)
    {
        grid_type mesh
            = fs::unstructured_mesh(points, indptr, indices, convex_hull_indices, areas, {});

        EXPECT_EQ(mesh.size(), mesh_size);
        EXPECT_EQ(mesh.shape(), shape_type({ mesh_size }));
    }

    TEST_F(unstructured_mesh, status_at_nodes)
    {
        {
            SCOPED_TRACE("default boundary conditions (convex hull nodes = fixed value)");

            grid_type mesh
                = fs::unstructured_mesh(points, indptr, indices, convex_hull_indices, areas, {});

            auto actual = mesh.status_at_nodes();

            EXPECT_EQ(actual(0), fs::node_status::fixed_value_boundary);
            EXPECT_EQ(actual(1), fs::node_status::fixed_value_boundary);
            EXPECT_EQ(actual(2), fs::node_status::fixed_value_boundary);
            EXPECT_EQ(actual(3), fs::node_status::fixed_value_boundary);
            EXPECT_EQ(actual(4), fs::node_status::core);
        }

        {
            SCOPED_TRACE("custom boundary conditions");

            grid_type mesh
                = fs::unstructured_mesh(points,
                                        indptr,
                                        indices,
                                        convex_hull_indices,
                                        areas,
                                        { { 2, fs::node_status::fixed_value_boundary } });

            auto actual = mesh.status_at_nodes();

            EXPECT_EQ(actual(0), fs::node_status::core);
            EXPECT_EQ(actual(1), fs::node_status::core);
            EXPECT_EQ(actual(2), fs::node_status::fixed_value_boundary);
            EXPECT_EQ(actual(3), fs::node_status::core);
            EXPECT_EQ(actual(4), fs::node_status::core);
        }

        {
            SCOPED_TRACE("looped boundary conditions not supported");

            EXPECT_THROW(fs::unstructured_mesh(points,
                                               indptr,
                                               indices,
                                               convex_hull_indices,
                                               areas,
                                               { { 2, fs::node_status::looped_boundary } }),
                         std::invalid_argument);
        }
    }

    TEST_F(unstructured_mesh, neighbors_count)
    {
        grid_type mesh
            = fs::unstructured_mesh(points, indptr, indices, convex_hull_indices, areas, {});

        EXPECT_EQ(mesh.neighbors_count(4), 3);  // 4
        EXPECT_EQ(mesh.neighbors_count(0), 3);
    }


    // test line 356 in base.hpp
    TEST_F(unstructured_mesh, neighbors_indices)  // only returns the indices - just the integer pos
    {
        grid_type mesh
            = fs::unstructured_mesh(points, indptr, indices, convex_hull_indices, areas, {});

        EXPECT_EQ(xt::all(xt::equal((mesh.neighbors(3)), 3)));
    }

    TEST_F(unstructured_mesh, neighbors_distances)
    {
        grid_type mesh
            = fs::unstructured_mesh(points, indptr, indices, convex_hull_indices, areas, {});


        EXPECT_EQ()
    }


    TEST_F(unstructured_mesh, neighbor)  // neighbor returns a struct with distance or status chech
                                         // base.hpp, the structs in base
    {
        grid_type mesh
            = fs::unstructured_mesh(points, indptr, indices, convex_hull_indices, areas, {});

        // index
        EXPECT_EQ(xt::all(xt::equal((mesh.neighbors(3)), 3)));
        // distance

        // status
    }
}
