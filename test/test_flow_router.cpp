#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
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
        class dummy_flow_router final : public fs::flow_router<FG>
        {
        public:
            using base_type = fs::flow_router<FG>;
            using elevation_type = typename base_type::elevation_type;

            dummy_flow_router() = default;

            virtual ~dummy_flow_router() = default;

            void route1(const elevation_type& /*elevation*/, FG& /*fgraph*/){};
            void route2(const elevation_type& /*elevation*/, FG& /*fgraph*/){};
        };


        class flow_router : public ::testing::Test
        {
        protected:
            using flow_graph_type = fs::flow_graph<fs::raster_grid>;
            using grid_type = fs::raster_grid;
            using size_type = typename grid_type::size_type;

            fs::node_status fb = fs::node_status::fixed_value_boundary;
            fs::raster_boundary_status fixed_value_status{ fb };

            grid_type grid = grid_type({ 4, 4 }, { 1.1, 1.2 }, fixed_value_status);

            xt::xarray<double> elevation{ { 0.82, 0.16, 0.14, 0.20 },
                                          { 0.71, 0.97, 0.41, 0.09 },
                                          { 0.49, 0.01, 0.19, 0.38 },
                                          { 0.29, 0.82, 0.09, 0.88 } };

            flow_graph_type graph
                = flow_graph_type(grid,
                                  std::make_unique<dummy_flow_router<flow_graph_type>>(),
                                  std::make_unique<fs::no_sink_resolver<flow_graph_type>>());

            void update()
            {
                graph.update_routes(elevation);
            };
        };


        TEST_F(flow_router, receivers)
        {
            update();

            EXPECT_TRUE(xt::all(xt::equal(graph.receivers(), xt::ones<int>({ 16, 8 }) * -1)));
        }


        TEST_F(flow_router, receivers_distance)
        {
            update();

            EXPECT_TRUE(xt::allclose(graph.receivers_distance(), xt::ones<double>({ 16, 8 }) * -1));
        }
    }
}
