/**
 * @brief Class used to resolve local depressions of
 * the topographic surface.
 *
 * ``sink_resolver`` is meant to be used in combination
 * with ``flow_router`` through a ``flow_graph``.
 *
 */
#ifndef FASTSCAPELIB_FLOW_SINK_RESOLVER_H
#define FASTSCAPELIB_FLOW_SINK_RESOLVER_H

#include "fastscapelib/flow/flow_graph_impl.hpp"


namespace fastscapelib
{

    namespace detail
    {

        /**
         * Common implementation for all sink resolving methods.
         *
         * @tparam FG The flow graph implementation type.
         * @tparam SR The sink resolver type.
         */
        template <class FG, class SR>
        class sink_resolver_impl_base
        {
        public:
            using graph_impl_type = FG;
            using resolver_type = SR;

            sink_resolver_impl_base(graph_impl_type& graph_impl, const resolver_type& resolver)
                : m_graph_impl(graph_impl)
                , m_resolver(resolver){};

        private:
            graph_impl_type& m_graph_impl;
            const resolver_type& m_resolver;
        };


        /**
         * Sink resolving implementation.
         *
         * This class is used via template specialization (one for each
         * flow routing method).
         *
         * The declaration for the generic case here contains the minimum that
         * should be (re)implemented in specialized template classes.
         *
         * @tparam FG The flow graph implementation type.
         * @tparam SR The sink resolver type.
         */
        template <class FG, class SR>
        class sink_resolver_impl : public sink_resolver_impl_base<FG, SR>
        {
        public:
            using graph_impl_type = FG;
            using resolver_type = SR;
            using base_type = sink_resolver_impl_base<graph_impl_type, resolver_type>;

            using elevation_type = typename graph_impl_type::elevation_type;

            sink_resolver_impl(graph_impl_type& graph, const resolver_type& resolver)
                : base_type(graph, resolver){};

            const elevation_type& resolve1(const elevation_type& elevation)
            {
                return elevation;
            };
            const elevation_type& resolve2(const elevation_type& elevation)
            {
                return elevation;
            };
        };
    }

    struct no_sink_resolver
    {
        using flow_graph_impl_tag = detail::flow_graph_fixed_array_tag;
    };
}

#endif
