#include <stdexcept>

#include "gtest/gtest.h"
#include "xtensor/xsort.hpp"
#include "xtensor/xindex_view.hpp"
#include "xtensor/xtensor.hpp"

#include "fastscapelib/grid/trimesh.hpp"


namespace fs = fastscapelib;


namespace fastscapelib
{
    namespace testing
    {

        class trimesh : public ::testing::Test
        {
        protected:
            using node_s = fs::node_status;
            node_s fixed = node_s::fixed_value;
            node_s core = node_s::core;

            using grid_type = fs::trimesh_xt<fs::xt_selector>;
            using size_type = typename grid_type::size_type;
            using shape_type = typename grid_type::shape_type;

            size_type mesh_size = 5;

            // See Python tests for details about the mesh arrays used here
            xt::xtensor<double, 2> points{
                { 0.0, 0.5 }, { 0.5, 0.0 }, { 0.0, -0.5 }, { -0.5, 0.0 }, { 0.0, 0.0 }
            };
            xt::xtensor<size_type, 2> triangles{
                { 2, 4, 3 }, { 4, 2, 1 }, { 4, 0, 3 }, { 0, 4, 1 }
            };
        };

        TEST_F(trimesh, static_expr)
        {
            EXPECT_EQ(fs::trimesh::is_structured(), false);
            EXPECT_EQ(fs::trimesh::is_uniform(), false);
            EXPECT_EQ(fs::trimesh::n_neighbors_max(), 20u);
            EXPECT_EQ(fs::trimesh::xt_ndims(), 1);
        }


        TEST_F(trimesh, ctor)
        {
            grid_type mesh = fs::trimesh(points, triangles);

            EXPECT_EQ(mesh.size(), mesh_size);
            EXPECT_EQ(mesh.shape(), shape_type({ mesh_size }));

            xt::xtensor<double, 2> invalid_points{ { 0., 1., 2. } };
            EXPECT_THROW(fs::trimesh(invalid_points, triangles), std::invalid_argument);

            xt::xtensor<double, 2> invalid_triangles{ { 0., 1. } };
            EXPECT_THROW(fs::trimesh(points, invalid_triangles), std::invalid_argument);

            xt::xtensor<node_status, 1> invalid_nodes_status{ fs::node_status::core };
            EXPECT_THROW(fs::trimesh(points, triangles, invalid_nodes_status),
                         std::invalid_argument);
        }

        TEST_F(trimesh, nodes_status)
        {
            {
                SCOPED_TRACE("default boundary conditions (boundary nodes = fixed value)");

                grid_type mesh = fs::trimesh(points, triangles);

                auto actual = mesh.nodes_status();

                EXPECT_EQ(actual(0), fs::node_status::fixed_value);
                EXPECT_EQ(actual(1), fs::node_status::fixed_value);
                EXPECT_EQ(actual(2), fs::node_status::fixed_value);
                EXPECT_EQ(actual(3), fs::node_status::fixed_value);
                EXPECT_EQ(actual(4), fs::node_status::core);
            }

            {
                SCOPED_TRACE("custom boundary conditions");

                grid_type mesh
                    = fs::trimesh(points, triangles, { { 2, fs::node_status::fixed_value } });
                auto actual = mesh.nodes_status();

                EXPECT_EQ(actual(0), fs::node_status::core);
                EXPECT_EQ(actual(1), fs::node_status::core);
                EXPECT_EQ(actual(2), fs::node_status::fixed_value);
                EXPECT_EQ(actual(3), fs::node_status::core);
                EXPECT_EQ(actual(4), fs::node_status::core);
            }

            {
                SCOPED_TRACE("custom boundary conditions (array)");

                auto core = fs::node_status::core;
                auto fixed = fs::node_status::fixed_value;

                xt::xtensor<node_status, 1> nodes_status{ core, core, fixed, core, core };
                grid_type mesh = fs::trimesh(points, triangles, nodes_status);

                EXPECT_EQ(mesh.nodes_status(), nodes_status);
            }

            {
                SCOPED_TRACE("looped boundary conditions not supported");

                EXPECT_THROW(fs::trimesh(points, triangles, { { 2, fs::node_status::looped } }),
                             std::invalid_argument);
            }
        }

        TEST_F(trimesh, nodes_areas)
        {
            // simple case of 1x1 square domain with 9 evenly spaced nodes
            xt::xtensor<double, 2> points{ { 0.0, 0.0 }, { 0.0, 0.5 }, { 0.0, 1.0 },
                                           { 0.5, 1.0 }, { 1.0, 1.0 }, { 1.0, 0.5 },
                                           { 1.0, 0.0 }, { 0.5, 0.0 }, { 0.5, 0.5 } };
            xt::xtensor<size_type, 2> triangles{
                { 0, 1, 8 }, { 1, 2, 8 }, { 2, 3, 8 }, { 3, 4, 8 },
                { 4, 5, 8 }, { 5, 6, 8 }, { 6, 7, 8 }, { 7, 0, 8 },
            };
            grid_type mesh = fs::trimesh(points, triangles);

            EXPECT_EQ(mesh.nodes_areas(0), 0.0625);
            EXPECT_EQ(mesh.nodes_areas(8), 0.25);

            EXPECT_EQ(mesh.nodes_areas(),
                      (xt::xtensor<double, 1>{
                          0.0625, 0.125, 0.0625, 0.125, 0.0625, 0.125, 0.0625, 0.125, 0.25 }));
        }

        TEST_F(trimesh, neighbors_count)
        {
            grid_type mesh = fs::trimesh(points, triangles);

            EXPECT_EQ(mesh.neighbors_count(0), 3);
            EXPECT_EQ(mesh.neighbors_count(4), 4);
        }

        TEST_F(trimesh, neighbors_indices)
        {
            grid_type mesh = fs::trimesh(points, triangles);

            EXPECT_EQ(xt::sort(mesh.neighbors_indices(3)),
                      (xt::xtensor<std::size_t, 1>{ 0, 2, 4 }));
            EXPECT_EQ(xt::sort(mesh.neighbors_indices(4)),
                      (xt::xtensor<std::size_t, 1>{ 0, 1, 2, 3 }));
        }

        TEST_F(trimesh, neighbors_distances)
        {
            grid_type mesh = fs::trimesh(points, triangles);

            double diag_dist = std::sqrt(0.5);

            {
                SCOPED_TRACE("node 3");
                auto sorted_nidx = xt::argsort(mesh.neighbors_indices(3));
                auto actual = xt::index_view(mesh.neighbors_distances(3), sorted_nidx);
                xt::xtensor<double, 1> expected{ diag_dist, diag_dist, 0.5 };

                EXPECT_TRUE(xt::allclose(actual, expected));
            }
            {
                SCOPED_TRACE("node 4");
                EXPECT_EQ(mesh.neighbors_distances(4),
                          (xt::xtensor<double, 1>{ 0.5, 0.5, 0.5, 0.5 }));
            }
        }

        TEST_F(trimesh, neighbor)
        {
            using neighbors_type = std::vector<fs::neighbor>;

            grid_type mesh = fs::trimesh(points, triangles);

            double diag_dist = std::sqrt(0.5);

            {
                SCOPED_TRACE("node 3");

                auto sorted_nidx = xt::argsort(mesh.neighbors_indices(3));
                auto temp = mesh.neighbors(3);
                EXPECT_EQ(temp.size(), sorted_nidx.size());

                neighbors_type actual(sorted_nidx.size());
                for (std::size_t i = 0; i < temp.size(); i++)
                {
                    actual[i] = temp[sorted_nidx(i)];
                }

                neighbors_type expected{ { 0, diag_dist, fixed },
                                         { 2, diag_dist, fixed },
                                         { 4, 0.5, core } };

                EXPECT_EQ(actual, expected);
            }
            {
                SCOPED_TRACE("node 4");

                auto sorted_nidx = xt::argsort(mesh.neighbors_indices(4));
                auto temp = mesh.neighbors(4);
                EXPECT_EQ(temp.size(), sorted_nidx.size());

                neighbors_type actual(sorted_nidx.size());
                for (std::size_t i = 0; i < temp.size(); i++)
                {
                    actual[i] = temp[sorted_nidx(i)];
                }

                neighbors_type expected{
                    { 0, 0.5, fixed }, { 1, 0.5, fixed }, { 2, 0.5, fixed }, { 3, 0.5, fixed }
                };

                EXPECT_EQ(actual, expected);
            }
        }
    }
}
