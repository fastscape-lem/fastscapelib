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


namespace fastscapelib
{

    struct sink_resolver
    {
        // enable dynamic downcasting on sink resolver types
        virtual ~sink_resolver() = default;
    };

    namespace detail
    {

        /**
         * Common implementation for all sink resolving methods.
         *
         * @tparam FG The flow graph type.
         * @tparam SR The sink resolver type.
         */
        template <class FG, class SR>
        class sink_resolver_impl_base
        {
        public:
            using graph_type = FG;
            using resolver_type = SR;

            sink_resolver_impl_base(graph_type& graph, const resolver_type& resolver)
                : m_graph(graph)
                , m_resolver(resolver){};

        private:
            graph_type& m_graph;
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
         * @tparam FG The flow graph type.
         * @tparam SR The sink resolver type.
         */
        template <class FG, class SR>
        class sink_resolver_impl : public sink_resolver_impl_base<FG, SR>
        {
        public:
            using graph_type = FG;
            using resolver_type = SR;
            using base_type = sink_resolver_impl_base<graph_type, resolver_type>;

            using elevation_type = typename graph_type::elevation_type;

            sink_resolver_impl(graph_type& graph, const resolver_type& resolver)
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

    struct no_sink_resolver : public sink_resolver
    {
    };
}

#endif
