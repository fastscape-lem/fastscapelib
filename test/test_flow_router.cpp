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

    struct test_flow_router
    {
        using flow_graph_impl_tag = detail::flow_graph_fixed_array_tag;
    };


    namespace detail
    {

        template <class FG>
        class flow_router_impl<FG, test_flow_router>
            : public flow_router_impl_base<FG, test_flow_router>
        {
        public:
            using graph_type = FG;
            using base_type = flow_router_impl_base<graph_type, test_flow_router>;

            using elevation_type = typename graph_type::elevation_type;

            static constexpr size_t n_receivers = 0;

            flow_router_impl(graph_type& graph, const test_flow_router& router)
                : base_type(graph, router){};

            void route1(const elevation_type& /*elevation*/){};
            void route2(const elevation_type& /*elevation*/){};
        };
    }

    namespace testing
    {

        class flow_router : public ::testing::Test
        {
        protected:
            using flow_graph_type
                = fs::flow_graph<fs::raster_grid, test_flow_router, fs::no_sink_resolver>;
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
                = fs::make_flow_graph(grid, fs::test_flow_router(), fs::no_sink_resolver());

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
