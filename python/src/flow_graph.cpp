#include <memory>

#include "xtensor-python/pytensor.hpp"

#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/stl.h"

#include "grid.hpp"
#include "flow_graph.hpp"

#include "fastscapelib/flow/flow_operator.hpp"
#include "fastscapelib/flow/flow_snapshot.hpp"

namespace py = pybind11;
namespace fs = fastscapelib;


namespace fastscapelib
{
    class op1 : public flow_operator
    {
    public:
        inline std::string name() const noexcept override
        {
            return "op1";
        }

        static constexpr bool elevation_updated = true;

        int param = 1;
    };

    class op2 : public flow_operator
    {
    public:
        inline std::string name() const noexcept override
        {
            return "op2";
        }

        static constexpr bool graph_updated = true;
        static const flow_direction out_flowdir = flow_direction::single;

        int param = 2;
    };

    namespace detail
    {
        template <class FG>
        class flow_operator_impl<FG, op1, fs::flow_graph_fixed_array_tag>
            : public flow_operator_impl_base<FG, op1>
        {
        public:
            using base_type = flow_operator_impl_base<FG, op1>;

            flow_operator_impl(std::shared_ptr<op1> ptr)
                : base_type(std::move(ptr)){};
        };

        template <class FG>
        class flow_operator_impl<FG, op2, fs::flow_graph_fixed_array_tag>
            : public flow_operator_impl_base<FG, op2>
        {
        public:
            using base_type = flow_operator_impl_base<FG, op2>;

            flow_operator_impl(std::shared_ptr<op2> ptr)
                : base_type(std::move(ptr)){};
        };
    }
}


void
add_flow_graph_bindings(py::module& m)
{
    /*
     * Flow graph implementation
     */

    py::class_<fs::py_flow_graph_impl>(m, "FlowGraphImpl")
        .def_property_readonly("receivers", &fs::py_flow_graph_impl::receivers)
        .def_property_readonly("receivers_count", &fs::py_flow_graph_impl::receivers_count)
        .def_property_readonly("receivers_distance", &fs::py_flow_graph_impl::receivers_distance)
        .def_property_readonly("receivers_weight", &fs::py_flow_graph_impl::receivers_weight)
        .def_property_readonly("donors", &fs::py_flow_graph_impl::donors)
        .def_property_readonly("donors_count", &fs::py_flow_graph_impl::donors_count)
        .def_property_readonly("dfs_indices", &fs::py_flow_graph_impl::dfs_indices)
        .def_property_readonly("basins", &fs::py_flow_graph_impl::basins);

    /*
     * Flow operators
     */

    py::class_<fs::flow_operator, std::shared_ptr<fs::flow_operator>>(m, "FlowOperator");

    py::class_<fs::op1, fs::flow_operator, std::shared_ptr<fs::op1>>(m, "Op1")
        .def(py::init<>())
        .def_readwrite("param", &fs::op1::param);
    py::class_<fs::op2, fs::flow_operator, std::shared_ptr<fs::op2>>(m, "Op2")
        .def(py::init<>())
        .def_readwrite("param", &fs::op2::param);

    py::class_<fs::flow_snapshot, fs::flow_operator, std::shared_ptr<fs::flow_snapshot>>(
        m, "FlowSnapshot")
        .def(py::init<std::string, bool, bool>(),
             py::arg("snapshot_name"),
             py::arg("save_graph") = true,
             py::arg("save_elevation") = false);

    /*
     * Flow graph
     */

    py::class_<fs::py_flow_graph> pyfgraph(m, "FlowGraph");

    fs::register_py_flow_graph_init<fs::py_profile_grid>(pyfgraph);
    fs::register_py_flow_graph_init<fs::py_raster_grid>(pyfgraph);

    pyfgraph.def("impl", &fs::py_flow_graph::impl, py::return_value_policy::reference);
    pyfgraph.def("update_routes", &fs::py_flow_graph::update_routes);

    using data_array_type = fs::py_flow_graph::data_array_type;
    using data_type = fs::py_flow_graph::data_type;

    pyfgraph
        .def("accumulate",
             py::overload_cast<data_array_type&, const data_array_type&>(
                 &fs::py_flow_graph::accumulate, py::const_))
        .def("accumulate",
             py::overload_cast<data_array_type&, data_type>(&fs::py_flow_graph::accumulate,
                                                            py::const_))
        .def("accumulate",
             py::overload_cast<const data_array_type&>(&fs::py_flow_graph::accumulate, py::const_))
        .def("accumulate",
             py::overload_cast<data_type>(&fs::py_flow_graph::accumulate, py::const_));

    pyfgraph.def("basins", &fs::py_flow_graph::basins);
}
