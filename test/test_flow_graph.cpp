#include "fastscapelib/flow_graph.hpp"
#include "fastscapelib/flow_router.hpp"
#include "fastscapelib/flow_router_factory.hpp"
#include "fastscapelib/sink_resolver.hpp"
#include "fastscapelib/raster_grid.hpp"

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
            
            constant_before_sink_resolver():
                m_elevation(elevation_type({0}))
            {}
            
            virtual ~constant_before_sink_resolver() = default;
            
            const elevation_type& resolve_before_route(const elevation_type& elevation, FG& /*fgraph*/) override
            { 
                m_elevation = xt::ones_like(elevation)*2.;
                return m_elevation;
            }
            
            const elevation_type& resolve_after_route(const elevation_type& elevation, FG& /*fgraph*/) override
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
            
            constant_after_sink_resolver():
                m_elevation(elevation_type({0}))
            {}
            
            virtual ~constant_after_sink_resolver() = default;
            
            const elevation_type& resolve_before_route(const elevation_type& elevation, FG& /*fgraph*/) override
            { 
                return elevation;
            }
            
            const elevation_type& resolve_after_route(const elevation_type& elevation, FG& /*fgraph*/) override
            {
                m_elevation = xt::ones_like(elevation)*5.;
                return m_elevation;
            }

        private:
            elevation_type m_elevation;
        };

       class flow_graph: public ::testing::Test
        {
        protected:

            using flow_graph_type = fs::flow_graph<fs::raster_grid, double>;
            using grid_type = fs::raster_grid;
            using size_type = typename grid_type::size_type;

            fs::node_status fb = fs::node_status::fixed_value_boundary;
            fs::raster_boundary_status fixed_value_status {fb};

            grid_type grid = grid_type({4, 4}, {1.1, 1.2}, fixed_value_status);

            xt::xtensor<double, 2> elevation
                {{0.82,  0.16,  0.14,  0.20},
                 {0.71,  0.97,  0.41,  0.09},
                 {0.49,  0.01,  0.19,  0.38},
                 {0.29,  0.82,  0.09,  0.88}};
        };

        TEST_F(flow_graph, ctor)
        {
            auto router = std::make_unique<fs::dummy_flow_router<flow_graph_type>>();
            auto resolver = std::make_unique<fs::no_sink_resolver<flow_graph_type>>();

            flow_graph_type graph(grid, std::move(router), std::move(resolver));

            EXPECT_EQ(graph.grid().size(), 16u); // dummy test
        }
    }
}