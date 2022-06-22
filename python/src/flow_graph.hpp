#ifndef PYFASTSCAPELIB_FLOW_GRAPH_H
#define PYFASTSCAPELIB_FLOW_GRAPH_H

#include <memory>

#include "pybind11/pybind11.h"

#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/flow/sink_resolver.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"

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
         * Flow graph facade class for Python bindings.
         *
         * It implements type erasure in order to expose
         * a single class to Python for all grid types.
         *
         */
        class py_flow_graph;


        class flow_graph_wrapper_base
        {
        public:
            using index_type = std::size_t;
            using neighbors_count_type = std::uint8_t;
            using grid_data_type = double;

            using data_type = xt_array_t<py_selector, double>;
            using donors_type = xt_tensor_t<py_selector, index_type, 2>;
            using donors_count_type = xt_tensor_t<py_selector, neighbors_count_type, 1>;

            using receivers_type = donors_type;
            using receivers_count_type = donors_count_type;
            using receivers_weight_type = xt_tensor_t<py_selector, double, 2>;
            using receivers_distance_type = xt_tensor_t<py_selector, grid_data_type, 2>;

            using stack_type = xt_tensor_t<py_selector, index_type, 1>;

            virtual ~flow_graph_wrapper_base(){};

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
            using flow_graph_type = fs::flow_graph<G, fs::py_selector>;

            using index_type = typename flow_graph_wrapper_base::index_type;
            using neighbors_count_type = typename flow_graph_wrapper_base::neighbors_count_type;
            using grid_data_type = typename flow_graph_wrapper_base::grid_data_type;
            using data_type = typename flow_graph_wrapper_base::data_type;
            using donors_type = typename flow_graph_wrapper_base::donors_type;
            using donors_count_type = typename flow_graph_wrapper_base::donors_count_type;
            using receivers_type = typename flow_graph_wrapper_base::receivers_type;
            using receivers_count_type = typename flow_graph_wrapper_base::receivers_count_type;
            using receivers_weight_type = typename flow_graph_wrapper_base::receivers_weight_type;
            using receivers_distance_type =
                typename flow_graph_wrapper_base::receivers_distance_type;
            using stack_type = typename flow_graph_wrapper_base::stack_type;

            flow_graph_wrapper(G& grid, py_flow_router& py_router, py_sink_resolver& py_resolver)
            {
                auto router_ptr = fs::detail::make_flow_router<flow_graph_type>(py_router);
                auto resolver_ptr = fs::detail::make_sink_resolver<flow_graph_type>(py_resolver);

                p_graph = std::make_unique<flow_graph_type>(
                    grid, std::move(router_ptr), std::move(resolver_ptr));
            }

            virtual ~flow_graph_wrapper(){};

            const data_type& update_routes(const data_type& elevation)
            {
                return p_graph->update_routes(elevation);
            };

            const receivers_type& receivers() const
            {
                return p_graph->receivers();
            };

            const receivers_count_type& receivers_count() const
            {
                return p_graph->receivers_count();
            };

            const receivers_distance_type& receivers_distance() const
            {
                return p_graph->receivers_distance();
            };

            const receivers_weight_type& receivers_weight() const
            {
                return p_graph->receivers_weight();
            };

            const donors_type& donors() const
            {
                return p_graph->donors();
            };

            const donors_count_type& donors_count() const
            {
                return p_graph->donors_count();
            };

            const stack_type& dfs_stack() const
            {
                return p_graph->dfs_stack();
            };

            data_type accumulate(const data_type& data) const
            {
                return p_graph->accumulate(data);
            };

            data_type accumulate(const double& data) const
            {
                return p_graph->accumulate(data);
            };

        private:
            std::unique_ptr<flow_graph_type> p_graph;
        };


        class py_flow_graph
        {
        public:
            using index_type = std::size_t;
            using neighbors_count_type = std::uint8_t;
            using grid_data_type = double;

            using data_type = xt_array_t<py_selector, double>;
            using donors_type = xt_tensor_t<py_selector, index_type, 2>;
            using donors_count_type = xt_tensor_t<py_selector, neighbors_count_type, 1>;

            using receivers_type = donors_type;
            using receivers_count_type = donors_count_type;
            using receivers_weight_type = xt_tensor_t<py_selector, double, 2>;
            using receivers_distance_type = xt_tensor_t<py_selector, grid_data_type, 2>;

            using stack_type = xt_tensor_t<py_selector, index_type, 1>;

            template <class G>
            py_flow_graph(G& grid, py_flow_router& py_router, py_sink_resolver& py_resolver)
                : p_wrapped_graph(
                    std::make_unique<flow_graph_wrapper<G>>(grid, py_router, py_resolver))
            {
            }

            const data_type& update_routes(const data_type& elevation)
            {
                return p_wrapped_graph->update_routes(elevation);
            };

            const receivers_type& receivers() const
            {
                return p_wrapped_graph->receivers();
            };

            const receivers_count_type& receivers_count() const
            {
                return p_wrapped_graph->receivers_count();
            };

            const receivers_distance_type& receivers_distance() const
            {
                return p_wrapped_graph->receivers_distance();
            };

            const receivers_weight_type& receivers_weight() const
            {
                return p_wrapped_graph->receivers_weight();
            };

            const donors_type& donors() const
            {
                return p_wrapped_graph->donors();
            };

            const donors_count_type& donors_count() const
            {
                return p_wrapped_graph->donors_count();
            };

            const stack_type& dfs_stack() const
            {
                return p_wrapped_graph->dfs_stack();
            };

            data_type accumulate(const data_type& data) const
            {
                return p_wrapped_graph->accumulate(data);
            };

            data_type accumulate(const double& data) const
            {
                return p_wrapped_graph->accumulate(data);
            };

        private:
            std::unique_ptr<flow_graph_wrapper_base> p_wrapped_graph;
        };
    }
}

#endif
