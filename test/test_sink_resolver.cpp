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

    struct test_sink_resolver
    {
    };

    namespace detail
    {

        template <class FG>
        class sink_resolver_impl<FG, test_sink_resolver>
            : public sink_resolver_impl_base<FG, test_sink_resolver>
        {
        public:
            using graph_type = FG;
            using base_type = sink_resolver_impl_base<graph_type, test_sink_resolver>;

            using elevation_type = typename graph_type::elevation_type;

            sink_resolver_impl(graph_type& graph, const test_sink_resolver& resolver)
                : base_type(graph, resolver)
                , m_elevation(elevation_type({ 0 })){};

            const elevation_type& resolve1(const elevation_type& elevation)
            {
                m_elevation = elevation + 10.;
                return m_elevation;
            };
            const elevation_type& resolve2(const elevation_type& elevation)
            {
                m_elevation = elevation + 5.;
                return m_elevation;
            };

        private:
            elevation_type m_elevation;
        };
    }

    namespace testing
    {

        class sink_resolver : public ::testing::Test
        {
        protected:
            using flow_graph_type
                = fs::flow_graph<fs::raster_grid, fs::single_flow_router, fs::test_sink_resolver>;
            using grid_type = fs::raster_grid;
            using size_type = typename grid_type::size_type;

            fs::node_status fb = fs::node_status::fixed_value_boundary;
            fs::raster_boundary_status fixed_value_status{ fb };

            grid_type grid = grid_type({ 4, 4 }, { 1.1, 1.2 }, fixed_value_status);

            xt::xtensor<double, 2> elevation{ { 0.82, 0.16, 0.14, 0.20 },
                                              { 0.71, 0.97, 0.41, 0.09 },
                                              { 0.49, 0.01, 0.19, 0.38 },
                                              { 0.29, 0.82, 0.09, 0.88 } };
        };


        TEST_F(sink_resolver, resolve)
        {
            auto graph
                = fs::make_flow_graph(grid, fs::single_flow_router(), fs::test_sink_resolver());

            const auto& new_elevation = graph.update_routes(elevation);
            EXPECT_TRUE(xt::all(xt::equal(new_elevation, elevation + 15.)));
        }
    }
}
