#ifndef FLOW_OPERATOR_FIXTURES_H_
#define FLOW_OPERATOR_FIXTURES_H_

#import "fastscapelib/flow/flow_operator.hpp"


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
}

#endif  // FLOW_OPERATOR_FIXTURES_H_
