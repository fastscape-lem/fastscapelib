#ifndef PYFASTSCAPELIB_FLOW_GRAPH_H
#define PYFASTSCAPELIB_FLOW_GRAPH_H

#include <memory>

#include "pybind11/pybind11.h"

#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_graph_impl.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/flow/sink_resolver.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"

#include "pytensor_utils.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


namespace fastscapelib
{

    /**
     * Flow graph implementation facade class for Python bindings.
     *
     * It implements type erasure in order to expose
     * a single class to Python for all grid types.
     *
     * It can be used to have read-only access to the graph underlying data
     * structures for any operation that is not supported by flow_graph.
     *
     * It is only accessed from a py_flow_graph object.
     *
     */
    class py_flow_graph_impl;


    namespace detail
    {

        class flow_graph_impl_wrapper_base
        {
        public:
            using neighbors_count_type = std::uint8_t;
            using grid_data_type = double;
            using index_type = std::size_t;

            using donors_type = xt_tensor_t<py_selector, index_type, 2>;
            using donors_count_type = xt_tensor_t<py_selector, neighbors_count_type, 1>;

            using receivers_type = donors_type;
            using receivers_count_type = donors_count_type;
            using receivers_weight_type = xt_tensor_t<py_selector, double, 2>;
            using receivers_distance_type = xt_tensor_t<py_selector, grid_data_type, 2>;

            using stack_type = xt_tensor_t<py_selector, index_type, 1>;

            virtual ~flow_graph_impl_wrapper_base(){};

            virtual const receivers_type& receivers() const = 0;

            virtual const receivers_count_type& receivers_count() const = 0;

            virtual const receivers_distance_type& receivers_distance() const = 0;

            virtual const receivers_weight_type& receivers_weight() const = 0;

            virtual const donors_type& donors() const = 0;

            virtual const donors_count_type& donors_count() const = 0;

            virtual const stack_type& dfs_stack() const = 0;
        };


        template <class FG>
        class flow_graph_impl_wrapper : public flow_graph_impl_wrapper_base
        {
        public:
            using flow_graph_impl_type = FG;

            using index_type = typename flow_graph_impl_wrapper_base::index_type;
            using neighbors_count_type =
                typename flow_graph_impl_wrapper_base::neighbors_count_type;
            using grid_data_type = typename flow_graph_impl_wrapper_base::grid_data_type;
            using donors_type = typename flow_graph_impl_wrapper_base::donors_type;
            using donors_count_type = typename flow_graph_impl_wrapper_base::donors_count_type;
            using receivers_type = typename flow_graph_impl_wrapper_base::receivers_type;
            using receivers_count_type =
                typename flow_graph_impl_wrapper_base::receivers_count_type;
            using receivers_weight_type =
                typename flow_graph_impl_wrapper_base::receivers_weight_type;
            using receivers_distance_type =
                typename flow_graph_impl_wrapper_base::receivers_distance_type;
            using stack_type = typename flow_graph_impl_wrapper_base::stack_type;

            virtual ~flow_graph_impl_wrapper(){};

            flow_graph_impl_wrapper(FG& graph_impl)
                : p_graph_impl(graph_impl){};

            const receivers_type& receivers() const
            {
                return p_graph_impl.receivers();
            };

            const receivers_count_type& receivers_count() const
            {
                return p_graph_impl.receivers_count();
            };

            const receivers_distance_type& receivers_distance() const
            {
                return p_graph_impl.receivers_distance();
            };

            const receivers_weight_type& receivers_weight() const
            {
                return p_graph_impl.receivers_weight();
            };

            const donors_type& donors() const
            {
                return p_graph_impl.donors();
            };

            const donors_count_type& donors_count() const
            {
                return p_graph_impl.donors_count();
            };

            const stack_type& dfs_stack() const
            {
                return p_graph_impl.dfs_stack();
            };

        private:
            flow_graph_impl_type& p_graph_impl;
        };
    }


    class py_flow_graph_impl
    {
    public:
        using index_type = std::size_t;
        using neighbors_count_type = std::uint8_t;
        using grid_data_type = double;

        using donors_type = xt_tensor_t<py_selector, index_type, 2>;
        using donors_count_type = xt_tensor_t<py_selector, neighbors_count_type, 1>;

        using receivers_type = donors_type;
        using receivers_count_type = donors_count_type;
        using receivers_weight_type = xt_tensor_t<py_selector, double, 2>;
        using receivers_distance_type = xt_tensor_t<py_selector, grid_data_type, 2>;

        using stack_type = xt_tensor_t<py_selector, index_type, 1>;


        template <class FG>
        py_flow_graph_impl(FG& graph_impl)
            : p_wrapped_graph_impl(
                std::make_unique<detail::flow_graph_impl_wrapper<FG>>(graph_impl)){};

        const receivers_type& receivers() const
        {
            return p_wrapped_graph_impl->receivers();
        };

        const receivers_count_type& receivers_count() const
        {
            return p_wrapped_graph_impl->receivers_count();
        };

        const receivers_distance_type& receivers_distance() const
        {
            return p_wrapped_graph_impl->receivers_distance();
        };

        const receivers_weight_type& receivers_weight() const
        {
            return p_wrapped_graph_impl->receivers_weight();
        };

        const donors_type& donors() const
        {
            return p_wrapped_graph_impl->donors();
        };

        const donors_count_type& donors_count() const
        {
            return p_wrapped_graph_impl->donors_count();
        };

        const stack_type& dfs_stack() const
        {
            return p_wrapped_graph_impl->dfs_stack();
        };

    private:
        std::unique_ptr<detail::flow_graph_impl_wrapper_base> p_wrapped_graph_impl;
    };


    /**
     * Flow graph facade class for Python bindings.
     *
     * It implements type erasure in order to expose
     * a single class to Python for all grid types.
     *
     */
    class py_flow_graph;


    namespace detail
    {

        class flow_graph_wrapper_base
        {
        public:
            using index_type = std::size_t;
            using data_type = xt_array_t<py_selector, double>;

            virtual ~flow_graph_wrapper_base(){};

            virtual index_type size() const = 0;

            virtual const py_flow_graph_impl& impl() const = 0;

            virtual const data_type& update_routes(const data_type& elevation) = 0;

            virtual data_type accumulate(const data_type& data) const = 0;

            virtual data_type accumulate(const double& data) const = 0;
        };

        template <class G, class FR, class SR>
        class flow_graph_wrapper : public flow_graph_wrapper_base
        {
        public:
            using flow_graph_type = fs::flow_graph<G, FR, SR, fs::py_selector>;

            using index_type = typename flow_graph_wrapper_base::index_type;
            using data_type = typename flow_graph_wrapper_base::data_type;

            flow_graph_wrapper(G& grid, const FR& router, const SR& resolver)
            {
                p_graph = std::make_unique<flow_graph_type>(grid, router, resolver);
                p_graph_impl = std::make_unique<py_flow_graph_impl>(p_graph->impl());
            }

            virtual ~flow_graph_wrapper(){};

            index_type size() const
            {
                return p_graph->size();
            };

            const py_flow_graph_impl& impl() const
            {
                return *p_graph_impl;
            };

            const data_type& update_routes(const data_type& elevation)
            {
                return p_graph->update_routes(elevation);
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
            std::unique_ptr<py_flow_graph_impl> p_graph_impl;
        };
    }


    class py_flow_graph
    {
    public:
        using index_type = std::size_t;
        using data_type = xt_array_t<py_selector, double>;

        template <class G, class FR, class SR>
        py_flow_graph(G& grid, const FR& router, const SR& resolver)
            : p_wrapped_graph(
                std::make_unique<detail::flow_graph_wrapper<G, FR, SR>>(grid, router, resolver)){};

        index_type size() const
        {
            return p_wrapped_graph->size();
        };

        const py_flow_graph_impl& impl() const
        {
            return p_wrapped_graph->impl();
        };

        const data_type& update_routes(const data_type& elevation)
        {
            return p_wrapped_graph->update_routes(elevation);
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
        std::unique_ptr<detail::flow_graph_wrapper_base> p_wrapped_graph;
    };


    template <class G>
    void add_init_methods(py::class_<py_flow_graph>& pyfg)
    {
        pyfg.def(py::init<G&, single_flow_router&, no_sink_resolver&>())
            .def(py::init<G&, multiple_flow_router&, no_sink_resolver&>());
    }
}

#endif
