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

    /**
     * Base class for the implementation of depression
     * filling or pit resolving.
     *
     * All derived classes must implement ``resolve1``
     * and ``resolve2`` methods.
     *
     * @tparam FG The flow_graph class.
     */
    template <class FG>
    class sink_resolver
    {
    public:
        using elevation_type = typename FG::elevation_type;

        // Entity semantics, i.e., a flow graph uses a sink resolver via this base class.
        // -> avoid incomplete destruction (e.g., there may be members in inherited classes
        // that need to be destroyed)
        virtual ~sink_resolver() = default;

        // do not assign or copy sink resolvers using the base class.
        sink_resolver(const sink_resolver&) = delete;
        sink_resolver(sink_resolver&&) = delete;
        sink_resolver& operator=(const sink_resolver&) = delete;
        sink_resolver& operator=(sink_resolver&&) = delete;

        virtual const elevation_type& resolve1(const elevation_type& elevation, FG& fgraph) = 0;
        virtual const elevation_type& resolve2(const elevation_type& elevation, FG& fgraph) = 0;

    protected:
        sink_resolver() = default;
    };

    /**
     * A sink resolver that doesn't resolve anything.
     *
     * @tparam FG The flow_graph class.
     */
    template <class FG>
    class no_sink_resolver final : public sink_resolver<FG>
    {
    public:
        using base_type = no_sink_resolver<FG>;
        using elevation_type = typename base_type::elevation_type;

        no_sink_resolver() = default;

        virtual ~no_sink_resolver() = default;

        const elevation_type& resolve1(const elevation_type& elevation, FG& /*fgraph*/) override
        {
            return elevation;
        }

        const elevation_type& resolve2(const elevation_type& elevation, FG& /*fgraph*/) override
        {
            return elevation;
        }
    };
}

#endif
