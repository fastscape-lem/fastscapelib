#ifndef PYFASTSCAPELIB_FLOW_GRAPH_H
#define PYFASTSCAPELIB_FLOW_GRAPH_H

#include <stdexcept>

#include "pybind11/pybind11.h"

#include "fastscapelib/flow_graph.hpp"
#include "fastscapelib/flow_router.hpp"
#include "fastscapelib/flow_router_factory.hpp"
#include "fastscapelib/sink_resolver.hpp"
#include "fastscapelib/sink_resolver_factory.hpp"
#include "fastscapelib/xtensor_utils.hpp"

#include "flow_router.hpp"
#include "sink_resolver.hpp"
#include "pytensor_utils.hpp"


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
            using index_type = std::size_t;
            using neighbors_count_type = std::uint8_t;
            using distance_type = double;

            using data_type = xt_container_t<pyarray_selector, double>;
            using donors_type = xt_container_t<pyarray_selector, index_type>;
            using donors_count_type = xt_container_t<pyarray_selector, neighbors_count_type>;

            using receivers_type = donors_type;
            using receivers_count_type = donors_count_type;
            using receivers_weight_type = xt_container_t<pyarray_selector, double>;
            using receivers_distance_type = xt_container_t<pyarray_selector, distance_type>;

            using stack_type = xt_container_t<pyarray_selector, index_type>;


            virtual ~flow_graph_wrapper_base(){};

            virtual int get_index_type() const = 0;

            virtual const data_type& update_routes(const data_type& elevation) = 0;

            virtual const receivers_type& receivers() const = 0;

            virtual const receivers_count_type& receivers_count() const = 0;

            virtual const receivers_distance_type& receivers_distance() const = 0;

            virtual const receivers_weight_type& receivers_weight() const = 0;

            virtual const donors_type& donors() const = 0;

            virtual const donors_count_type& donors_count() const = 0;

            virtual const stack_type& dfs_stack() const = 0;

            virtual data_type accumulate(const data_type& data) const = 0;

            virtual data_type accumulate(const double& data) const = 0;
        };

        template <class G>
        class flow_graph_wrapper : public flow_graph_wrapper_base
        {
        public:
            using wrapped_type = fs::flow_graph<G, double, fs::pyarray_selector>;

            using index_type = typename flow_graph_wrapper_base::index_type;
            using neighbors_count_type = typename flow_graph_wrapper_base::neighbors_count_type;
            using distance_type = typename flow_graph_wrapper_base::distance_type;
            using data_type = typename flow_graph_wrapper_base::data_type;
            using donors_type = typename flow_graph_wrapper_base::donors_type;
            using donors_count_type = typename flow_graph_wrapper_base::donors_count_type;
            using receivers_type = typename flow_graph_wrapper_base::receivers_type;
            using receivers_count_type = typename flow_graph_wrapper_base::receivers_count_type;
            using receivers_weight_type = typename flow_graph_wrapper_base::receivers_weight_type;
            using receivers_distance_type =
                typename flow_graph_wrapper_base::receivers_distance_type;
            using stack_type = typename flow_graph_wrapper_base::stack_type;

            flow_graph_wrapper(G& grid,
                               std::shared_ptr<flow_router_method> router,
                               std::shared_ptr<sink_resolver_method> resolver)
            {
                auto router_ptr = fs::detail::flow_router_factory<wrapped_type>::build(
                    router->method, *router->parameters);
                auto resolver_ptr
                    = fs::detail::sink_resolver_factory<wrapped_type>::build(resolver->method);

                if (router_ptr == nullptr)
                {
                    throw std::runtime_error("Using an unregistered flow router builder.");
                }
                if (resolver_ptr == nullptr)
                {
                    throw std::runtime_error("Using an unregistered sink resolver builder.");
                }

                m_wrapped = std::make_unique<wrapped_type>(
                    grid, std::move(router_ptr), std::move(resolver_ptr));
            }

            virtual ~flow_graph_wrapper(){};

            wrapped_type& get_wrapped()
            {
                return *m_wrapped;
            };

            int get_index_type() const
            {
                return 0;
            };

            const data_type& update_routes(const data_type& elevation)
            {
                return m_wrapped->update_routes(elevation);
            };

            const receivers_type& receivers() const
            {
                return m_wrapped->receivers();
            };

            const receivers_count_type& receivers_count() const
            {
                return m_wrapped->receivers_count();
            };

            const receivers_distance_type& receivers_distance() const
            {
                return m_wrapped->receivers_distance();
            };

            const receivers_weight_type& receivers_weight() const
            {
                return m_wrapped->receivers_weight();
            };

            const donors_type& donors() const
            {
                return m_wrapped->donors();
            };

            const donors_count_type& donors_count() const
            {
                return m_wrapped->donors_count();
            };

            const stack_type& dfs_stack() const
            {
                return m_wrapped->dfs_stack();
            };

            data_type accumulate(const data_type& data) const
            {
                return m_wrapped->accumulate(data);
            };

            data_type accumulate(const double& data) const
            {
                return m_wrapped->accumulate(data);
            };

        private:
            std::unique_ptr<wrapped_type> m_wrapped;
        };


        class flow_graph_facade
        {
        public:
            using self_type = flow_graph_facade;

            using index_type = std::size_t;
            using neighbors_count_type = std::uint8_t;
            using distance_type = double;

            using data_type = xt_container_t<pyarray_selector, double>;
            using donors_type = xt_container_t<pyarray_selector, index_type>;
            using donors_count_type = xt_container_t<pyarray_selector, neighbors_count_type>;

            using receivers_type = donors_type;
            using receivers_count_type = donors_count_type;
            using receivers_weight_type = xt_container_t<pyarray_selector, double>;
            using receivers_distance_type = xt_container_t<pyarray_selector, distance_type>;

            using stack_type = xt_container_t<pyarray_selector, index_type>;

            template <class G>
            flow_graph_facade(G& obj,
                              std::shared_ptr<flow_router_method> router,
                              std::shared_ptr<sink_resolver_method> resolver)
                : p_impl(std::make_unique<flow_graph_wrapper<G>>(obj, router, resolver))
            {
            }

            template <class G>
            fs::flow_graph<G, double, fs::pyarray_selector>& get_implementation()
            {
                auto& derived = dynamic_cast<flow_graph_wrapper<G>&>(*p_impl);
                return derived.get_wrapped();
            };

            int get_implementation_index_type() const
            {
                return p_impl->get_index_type();
            }

            const data_type& update_routes(const data_type& elevation)
            {
                return p_impl->update_routes(elevation);
            };

            const receivers_type& receivers() const
            {
                return p_impl->receivers();
            };

            const receivers_count_type& receivers_count() const
            {
                return p_impl->receivers_count();
            };

            const data_type& receivers_distance() const
            {
                return p_impl->receivers_distance();
            };

            const receivers_weight_type& receivers_weight() const
            {
                return p_impl->receivers_weight();
            };

            const donors_type& donors() const
            {
                return p_impl->donors();
            };

            const donors_count_type& donors_count() const
            {
                return p_impl->donors_count();
            };

            const stack_type& dfs_stack() const
            {
                return p_impl->dfs_stack();
            };

            data_type accumulate(const data_type& data) const
            {
                return p_impl->accumulate(data);
            };

            data_type accumulate(const double& data) const
            {
                return p_impl->accumulate(data);
            };

        private:
            std::unique_ptr<flow_graph_wrapper_base> p_impl;
        };
    }
}

#endif
