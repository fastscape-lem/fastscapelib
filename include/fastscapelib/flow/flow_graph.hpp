/**
 * @brief Class used to compute or follow flow routes on
 * the topographic surface.
 *
 * It provides a single API to chain flow_router and
 * sink_resolver methods and access their results.
 *
 */
#ifndef FASTSCAPELIB_FLOW_FLOW_GRAPH_H
#define FASTSCAPELIB_FLOW_FLOW_GRAPH_H

#include <stdexcept>
#include <type_traits>

#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/flow/sink_resolver.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"


namespace fastscapelib
{

    /**
     * Main class used to compute or follow flow routes on
     * the topographic surface.
     *
     * @tparam G The grid type.
     * @tparam FR The flow router type.
     * @tparam SR The sink resolver type.
     * @tparam S The xtensor selector type.
     */
    template <class G, class FR, class SR, class S = typename G::xt_selector>
    class flow_graph
    {
    public:
        using self_type = flow_graph<G, FR, SR, S>;
        using grid_type = G;
        using router_type = FR;
        using resolver_type = SR;
        using xt_selector = S;

        static_assert(std::is_same<typename router_type::flow_graph_impl_tag,
                                   typename resolver_type::flow_graph_impl_tag>::value,
                      "incompatible flow router and sink resolver types");
        using flow_graph_impl_tag = typename router_type::flow_graph_impl_tag;
        using flow_graph_impl_type
            = detail::flow_graph_impl<grid_type, xt_selector, flow_graph_impl_tag>;

        using size_type = typename grid_type::size_type;
        using elevation_type = xt_array_t<xt_selector, typename grid_type::grid_data_type>;

        template <class T>
        using data_type = xt_array_t<xt_selector, T>;

        flow_graph(G& grid, const router_type& router, const resolver_type& resolver)
            : m_grid(grid)
            , m_graph_impl(grid)
            , m_router_impl(m_graph_impl, router)
            , m_resolver_impl(m_graph_impl, resolver){};

        const elevation_type& update_routes(const elevation_type& elevation)
        {
            const auto& modified_elevation = m_resolver_impl.resolve1(elevation);
            m_router_impl.route1(modified_elevation);
            const auto& final_elevation = m_resolver_impl.resolve2(modified_elevation);
            m_router_impl.route2(final_elevation);

            return final_elevation;
        }

        G& grid()
        {
            return m_grid;
        };

        size_type size() const
        {
            return m_grid.size();
        };

        const flow_graph_impl_type& impl() const
        {
            return m_graph_impl;
        }

        template <class T>
        T accumulate(const T& data) const;

        data_type<double> accumulate(const double& data) const;

    private:
        G& m_grid;

        flow_graph_impl_type m_graph_impl;

        using router_impl_type =
            typename detail::flow_router_impl<flow_graph_impl_type, router_type>;
        using resolver_impl_type =
            typename detail::sink_resolver_impl<flow_graph_impl_type, resolver_type>;

        router_impl_type m_router_impl;
        resolver_impl_type m_resolver_impl;
    };

    template <class G, class FR, class SR, class S>
    template <class T>
    auto flow_graph<G, FR, SR, S>::accumulate(const T& data) const -> T
    {
        return m_graph_impl.accumulate(data);
    }

    template <class G, class FR, class SR, class S>
    auto flow_graph<G, FR, SR, S>::accumulate(const double& data) const -> data_type<double>
    {
        return m_graph_impl.accumulate(data);
    }


    template <class G, class FR, class SR, class S = typename G::xt_selector>
    flow_graph<G, FR, SR, S> make_flow_graph(G& grid, const FR& router, const SR& resolver)
    {
        return flow_graph<G, FR, SR, S>(grid, router, resolver);
    }
}

#endif
