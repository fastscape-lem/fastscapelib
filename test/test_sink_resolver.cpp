#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/flow/flow_router_factory.hpp"
#include "fastscapelib/flow/sink_resolver.hpp"
#include "fastscapelib/grid/raster_grid.hpp"

#include "xtensor/xtensor.hpp"
#include "xtensor/xio.hpp"

#include "gtest/gtest.h"


namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace testing
    {

        template <class FG>
        class constant_before_sink_resolver final : public fs::sink_resolver<FG>
        {
        public:
            using base_type = constant_before_sink_resolver<FG>;
            using elevation_type = typename base_type::elevation_type;

            constant_before_sink_resolver()
                : m_elevation(elevation_type({ 0 }))
            {
            }

            virtual ~constant_before_sink_resolver() = default;

            const elevation_type& resolve1(const elevation_type& elevation, FG& /*fgraph*/) override
            {
                m_elevation = xt::ones_like(elevation) * 2.;
                return m_elevation;
            }

            const elevation_type& resolve2(const elevation_type& elevation, FG& /*fgraph*/) override
            {
                return elevation;
            }

        private:
            elevation_type m_elevation;
        };

        template <class FG>
        class constant_after_sink_resolver final : public fs::sink_resolver<FG>
        {
        public:
            using base_type = constant_after_sink_resolver<FG>;
            using elevation_type = typename base_type::elevation_type;

            constant_after_sink_resolver()
                : m_elevation(elevation_type({ 0 }))
            {
            }

            virtual ~constant_after_sink_resolver() = default;

            const elevation_type& resolve1(const elevation_type& elevation, FG& /*fgraph*/) override
            {
                return elevation;
            }

            const elevation_type& resolve2(const elevation_type& elevation, FG& /*fgraph*/) override
            {
                m_elevation = xt::ones_like(elevation) * 5.;
                return m_elevation;
            }

        private:
            elevation_type m_elevation;
        };

        class sink_resolver : public ::testing::Test
        {
        protected:
            using flow_graph_type = fs::flow_graph<fs::raster_grid, double>;
            using grid_type = fs::raster_grid;
            using size_type = typename grid_type::size_type;

            fs::node_status fb = fs::node_status::fixed_value_boundary;
            fs::raster_boundary_status fixed_value_status{ fb };

            grid_type grid = grid_type({ 4, 4 }, { 1.1, 1.2 }, fixed_value_status);

            xt::xtensor<double, 2> elevation{ { 0.82, 0.16, 0.14, 0.20 },
                                              { 0.71, 0.97, 0.41, 0.09 },
                                              { 0.49, 0.01, 0.19, 0.38 },
                                              { 0.29, 0.82, 0.09, 0.88 } };

            xt::xtensor<size_type, 1> receivers = xt::ones<size_type>({ 16 }) * -1;
            xt::xtensor<size_type, 1> expected_receivers{ 1, 2, 7, 7, 9, 9, 7, 7,
                                                          9, 9, 9, 7, 9, 9, 9, 14 };

            double dia = std::sqrt(1.1 * 1.1 + 1.2 * 1.2);
            xt::xtensor<double, 1> dist2receivers = xt::ones<double>({ 16 }) * -1.;
            xt::xtensor<double, 1> expected_dist2receivers{
                1.2, 1.2, dia, 1.1, dia, 1.1, 1.2, 0.0, 1.2, 0.0, 1.2, 1.1, dia, 1.1, dia, 1.2
            };
        };


        TEST_F(sink_resolver, ctor)
        {
            fs::no_sink_resolver<flow_graph_type> sink();
        }


        TEST_F(sink_resolver, resolve1)
        {
            auto router = std::make_unique<fs::single_flow_router<flow_graph_type>>();
            auto resolver = std::make_unique<constant_before_sink_resolver<flow_graph_type>>();
            flow_graph_type graph(grid, std::move(router), std::move(resolver));

            xt::xarray<int> dummy_receivers = xt::ones<int>({ 16, 8 }) * -1;
            xt::view(dummy_receivers, xt::all(), 0) = xt::arange<int>(0, 16, 1);
            xt::xarray<double> dummy_dist2receivers = xt::ones<double>({ 16, 8 }) * -1;
            xt::view(dummy_dist2receivers, xt::all(), 0) = xt::zeros<double>({ 16 });

            const auto& graph_elevation = graph.update_routes(elevation);
            EXPECT_TRUE(xt::all(xt::equal(xt::ones_like(elevation) * 2., graph_elevation)));
            EXPECT_TRUE(xt::all(xt::equal(graph.receivers(), dummy_receivers)));
            EXPECT_TRUE(xt::allclose(graph.receivers_distance(), dummy_dist2receivers));
        }


        TEST_F(sink_resolver, resolve2)
        {
            auto router = std::make_unique<fs::single_flow_router<flow_graph_type>>();
            auto resolver = std::make_unique<constant_after_sink_resolver<flow_graph_type>>();
            auto graph = flow_graph_type(grid, std::move(router), std::move(resolver));

            const auto& graph_elevation = graph.update_routes(elevation);
            EXPECT_TRUE(xt::all(xt::equal(xt::ones_like(elevation) * 5., graph_elevation)));
            EXPECT_TRUE(
                xt::all(xt::equal(xt::view(graph.receivers(), xt::all(), 0), expected_receivers)));
            EXPECT_TRUE(xt::allclose(xt::view(graph.receivers_distance(), xt::all(), 0),
                                     expected_dist2receivers));
        }
    }
}
