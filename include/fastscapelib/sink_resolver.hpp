/**
 * @brief Class used to resolve local depressions of
 * the topographic surface.
 *
 * ``sink_resolver`` is meant to be used in combination
 * with ``flow_router`` through a ``flow_graph``.
 *
 */
#ifndef FASTSCAPELIB_SINK_RESOLVER_H
#define FASTSCAPELIB_SINK_RESOLVER_H


namespace fastscapelib
{

    /**
     * Base class for the implementation of depression
     * filling or pit resolving.
     *
     * All derived classes must implement ``resolve_before_route``
     * and ``resolve_after_route`` methods.
     *
     * @tparam FG The flow_graph class.
     */
    template <class FG>
    class sink_resolver
    {
    public:
        using elevation_type = typename FG::elevation_type;

        // Entity semantic
        virtual ~sink_resolver() = default;

        sink_resolver(const sink_resolver&) = delete;
        sink_resolver(sink_resolver&&) = delete;
        sink_resolver& operator=(const sink_resolver&) = delete;
        sink_resolver& operator=(sink_resolver&&) = delete;

        virtual const elevation_type& resolve_before_route(const elevation_type& elevation,
                                                           FG& fgraph)
            = 0;
        virtual const elevation_type& resolve_after_route(const elevation_type& elevation,
                                                          FG& fgraph)
            = 0;

    protected:
        sink_resolver() = default;
    };


    /**
     * A sink resolver with no action on the
     * the topographic surface.
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

        const elevation_type& resolve_before_route(const elevation_type& elevation,
                                                   FG& /*fgraph*/) override
        {
            return elevation;
        }

        const elevation_type& resolve_after_route(const elevation_type& elevation,
                                                  FG& /*fgraph*/) override
        {
            return elevation;
        }
    };


    /**
     * The possible sink resolvers.
     *
     */
    enum class sink_resolver_methods
    {
        none = 0,
        fill_pflood,
        fill_mst_kruskal,
        fill_mst_boruvka,
        fill_auto,  // either pflood or mst depending on number of sinks
        carve_mst_kruskal,
        carve_mst_boruvka
    };
}

#endif
