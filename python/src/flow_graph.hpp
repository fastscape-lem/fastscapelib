#ifndef PYFASTSCAPELIB_FLOW_GRAPH_H
#define PYFASTSCAPELIB_FLOW_GRAPH_H

#include "flow_router.hpp"
#include "sink_resolver.hpp"
#include "pytensor_utils.hpp"

#include "fastscapelib/flow_graph.hpp"
#include "fastscapelib/flow_router.hpp"
#include "fastscapelib/flow_router_factory.hpp"
#include "fastscapelib/sink_resolver.hpp"
#include "fastscapelib/sink_resolver_factory.hpp"

#include "pybind11/pybind11.h"

#include <stdexcept>


namespace py = pybind11;
namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace detail
    {
        /** 
         * Implementation of type erasure for the 
         * ``flow_graph`` class.
         * 
         * This allow to expose a single class through Python
         * bindings.
         * 
         */

        class flow_graph_facade;


        class flow_graph_wrapper_base
        {
        public:

            using elevation_type = xt::pyarray<double>;

            virtual ~flow_graph_wrapper_base() {};

            virtual int get_index_type() const = 0;

            virtual const elevation_type& update_routes(const elevation_type& elevation) = 0;
            
            virtual const xt::pyarray<std::size_t>& receivers() = 0;

            virtual const xt::pyarray<std::uint8_t>& receivers_count() = 0;

            virtual const xt::pyarray<double>& receivers_distance() = 0;

            virtual const xt::pyarray<double>& receivers_weight() = 0;
                        
            virtual const xt::pyarray<std::size_t>& donors() = 0;

            virtual const xt::pyarray<std::uint8_t>& donors_count() = 0;

            virtual const xt::pyarray<std::size_t>& dfs_stack() = 0;
        };

        template <class G>
        class flow_graph_wrapper : public flow_graph_wrapper_base
        {
        public:

            using wrapped_type = fs::flow_graph<G, double, fs::pyarray_selector>;
            using elevation_type = typename flow_graph_wrapper_base::elevation_type;

            flow_graph_wrapper(G& grid, std::shared_ptr<flow_router_method> router, std::shared_ptr<sink_resolver_method> resolver)
            {
                auto router_ptr = fs::detail::flow_router_factory<wrapped_type>::build(router->method, 
                                                                                       *router->parameters);
                auto resolver_ptr = fs::detail::sink_resolver_factory<wrapped_type>::build(resolver->method);

                if (router_ptr == nullptr) { throw std::runtime_error("Using an unregistered flow router builder."); }
                if (resolver_ptr == nullptr) { throw std::runtime_error("Using an unregistered sink resolver builder."); }

                m_wrapped = std::make_unique<wrapped_type>(grid,
                                                           std::move(router_ptr),
                                                           std::move(resolver_ptr));
            }

            virtual ~flow_graph_wrapper() {};

            wrapped_type& get_wrapped() { return *m_wrapped; };
            
            int get_index_type() const { return 0; };

            const elevation_type& update_routes(const elevation_type& elevation) { return m_wrapped->update_routes(elevation); };

            const xt::pyarray<std::size_t>& receivers() { return m_wrapped->receivers(); };

            const xt::pyarray<std::uint8_t>& receivers_count() { return m_wrapped->receivers_count(); };

            const xt::pyarray<double>& receivers_distance() { return m_wrapped->receivers_distance(); };

            const xt::pyarray<double>& receivers_weight() { return m_wrapped->receivers_weight(); };

            const xt::pyarray<std::size_t>& donors() { return m_wrapped->donors(); };

            const xt::pyarray<std::uint8_t>& donors_count() { return m_wrapped->donors_count(); };

            const xt::pyarray<std::size_t>& dfs_stack() { return m_wrapped->dfs_stack(); };

        private:

            std::unique_ptr<wrapped_type> m_wrapped;
        };


        class flow_graph_facade
        {
        public:

            using self_type = flow_graph_facade;
            using elevation_type = xt::pyarray<double>;

            template <class G>
            flow_graph_facade(G& obj, std::shared_ptr<flow_router_method> router, std::shared_ptr<sink_resolver_method> resolver)
                : p_impl(std::make_unique<flow_graph_wrapper<G>>(obj, router, resolver))
                {}

            template <class G>
            fs::flow_graph<G, elevation_type> get_implementation()
            {
                auto& derived = dynamic_cast<flow_graph_wrapper<G>&>(*p_impl);
                return derived.get_wrapped();
            };

            int get_implementation_index_type() const { return p_impl->get_index_type(); }

            const elevation_type& update_routes(const elevation_type& elevation) { return p_impl->update_routes(elevation); };

            const xt::pyarray<std::size_t>& receivers() { return p_impl->receivers(); };

            const xt::pyarray<std::uint8_t>& receivers_count() { return p_impl->receivers_count(); };

            const xt::pyarray<double>& receivers_distance() { return p_impl->receivers_distance(); };

            const xt::pyarray<double>& receivers_weight() { return p_impl->receivers_weight(); };

            const xt::pyarray<std::size_t>& donors() { return p_impl->donors(); };

            const xt::pyarray<std::uint8_t>& donors_count() { return p_impl->donors_count(); };

            const xt::pyarray<std::size_t>& dfs_stack() { return p_impl->dfs_stack(); };

        private:

            std::unique_ptr<flow_graph_wrapper_base> p_impl;
        };
    }
}

#endif