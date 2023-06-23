#include <memory>

#include "xtensor/xarray.hpp"
#include "xtensor/xbuilder.hpp"

#include "gtest/gtest.h"

#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_operator.hpp"
#include "fastscapelib/grid/raster_grid.hpp"


namespace fs = fastscapelib;

namespace fastscapelib
{

    /*
     * A flow operator for testing default values.
     */
    class default_operator : public flow_operator
    {
    public:
        inline std::string name() const noexcept override
        {
            return "default_operator";
        }
    };

    /*
     * A flow operator that does nothing (for testing).
     */
    class fake_operator : public flow_operator
    {
    public:
        inline std::string name() const noexcept override
        {
            return "fake_operator";
        }

        // does not update the graph but set to true so we can actually create a
        // flow graph without any error thrown.
        static constexpr bool graph_updated = true;
        static constexpr flow_direction out_flowdir = flow_direction::single;
    };

    namespace detail
    {

        template <class FG, class Tag>
        class flow_operator_impl<FG, fake_operator, Tag>
            : public flow_operator_impl_base<FG, fake_operator>
        {
        public:
            using base_type = flow_operator_impl_base<FG, fake_operator>;

            flow_operator_impl(std::shared_ptr<fake_operator> ptr)
                : base_type(std::move(ptr)){};
        };
    }

    /*
     * A trivial flow operator for testing.
     */
    class test_operator : public flow_operator
    {
    public:
        inline std::string name() const noexcept override
        {
            return "test_operator";
        }

        double m_offset = 1.0;

        static constexpr bool elevation_updated = true;
    };

    namespace detail
    {

        template <class FG, class Tag>
        class flow_operator_impl<FG, test_operator, Tag>
            : public flow_operator_impl_base<FG, test_operator>
        {
        public:
            using graph_impl_type = FG;
            using base_type = flow_operator_impl_base<graph_impl_type, test_operator>;

            using data_array_type = typename graph_impl_type::data_array_type;

            flow_operator_impl(std::shared_ptr<test_operator> ptr)
                : base_type(std::move(ptr)){};

            void apply(graph_impl_type& /*graph_impl*/, data_array_type& elevation)
            {
                elevation += this->m_op_ptr->m_offset;
            }
        };
    }

    namespace testing
    {

        TEST(flow_operator, default_member_values)
        {
            default_operator op;

            EXPECT_EQ(op.name(), "default_operator");
            EXPECT_EQ(op.graph_updated, false);
            EXPECT_EQ(op.elevation_updated, false);
            EXPECT_EQ(op.in_flowdir, fs::flow_direction::undefined);
            EXPECT_EQ(op.out_flowdir, fs::flow_direction::undefined);
        }

        TEST(flow_operator, return_same_elevation_ref)
        {
            fs::raster_boundary_status bs{ fs::node_status::fixed_value };
            auto grid = fs::raster_grid({ 5, 5 }, { 1.0, 1.0 }, bs);
            xt::xarray<double> elevation = xt::zeros<double>(grid.shape());

            test_operator op;
            auto graph = fs::flow_graph<fs::raster_grid>(grid, { fake_operator() });

            const auto& actual = graph.update_routes(elevation);

            // when elevation is not updated by any of the operators, update_routes
            // should return a const reference of the input elevation.
            ASSERT_TRUE(&actual == &elevation);
        }

        TEST(flow_operator, apply)
        {
            fs::raster_boundary_status bs{ fs::node_status::fixed_value };
            auto grid = fs::raster_grid({ 5, 5 }, { 1.0, 1.0 }, bs);
            xt::xarray<double> elevation = xt::zeros<double>(grid.shape());

            auto op = std::make_shared<test_operator>();
            auto graph = fs::flow_graph<fs::raster_grid>(grid, { op, fake_operator() });

            {
                SCOPED_TRACE("test apply 1st pass");

                const auto& actual = graph.update_routes(elevation);
                xt::xarray<double> expected = elevation + op->m_offset;

                EXPECT_EQ(actual, expected);
            }

            {
                SCOPED_TRACE("test apply 2nd pass (update operator params)");

                op->m_offset = 2.0;
                const auto& actual = graph.update_routes(elevation);
                xt::xarray<double> expected = elevation + op->m_offset;

                EXPECT_EQ(actual, expected);
            }
        }
    }
}
