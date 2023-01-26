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


template <class Tag>
struct flow_graph_impl_test
{
    using impl_tag = Tag;
    using data_array_type = xt::pyarray<double>;
};

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
        class flow_operator_impl<FG, op1, fs::detail::flow_graph_fixed_array_tag>
            : public flow_operator_impl_base<FG, op1>
        {
        public:
            using base_type = flow_operator_impl_base<FG, op1>;

            flow_operator_impl(std::shared_ptr<op1> ptr)
                : base_type(std::move(ptr)){};
        };

        template <class FG>
        class flow_operator_impl<FG, op2, fs::detail::flow_graph_fixed_array_tag>
            : public flow_operator_impl_base<FG, op2>
        {
        public:
            using base_type = flow_operator_impl_base<FG, op2>;

            flow_operator_impl(std::shared_ptr<op2> ptr)
                : base_type(std::move(ptr)){};
        };
    }

    /*
     * Only for Python bindings (need to add operators one at time from a py::list object).
     */
    template <class FG, class OPs>
    flow_operator_sequence<FG> make_flow_operator_sequence(OPs&& ops)
    {
        flow_operator_sequence<FG> op_sequence;

        for (auto op : ops)
        {
            // TODO: smarter way to do this? fold expressions?
            try
            {
                op_sequence.add_operator(op.template cast<std::shared_ptr<op1>>());
                continue;
            }
            catch (py::cast_error e)
            {
            }
            try
            {
                op_sequence.add_operator(op.template cast<std::shared_ptr<op2>>());
                continue;
            }
            catch (py::cast_error e)
            {
            }
            try
            {
                op_sequence.add_operator(op.template cast<std::shared_ptr<flow_snapshot>>());
                continue;
            }
            catch (py::cast_error e)
            {
                throw py::type_error("invalid flow operator");
            }
        }

        return std::move(op_sequence);
    }
}


class flow_graph_test
{
public:
    using impl_type = flow_graph_impl_test<fs::detail::flow_graph_fixed_array_tag>;
    using data_array_type = impl_type::data_array_type;

    using graph_impl_map = std::map<std::string, impl_type>;
    using elevation_map = std::map<std::string, data_array_type>;

    flow_graph_test() = default;
    flow_graph_test(fs::flow_operator_sequence<impl_type> operators)
        : m_operators(std::move(operators))
    {
        if (!m_operators.graph_updated())
        {
            throw std::invalid_argument(
                "must have at least one operator that updates the flow graph");
        }
        if (m_operators.out_flowdir() == fs::flow_direction::undefined)
        {
            throw std::invalid_argument(
                "must have at least one operator that defines the flow direction type");
        }

        // pre-allocate graph and elevation snapshots
        for (const auto& key : operators.graph_snapshot_keys())
        {
            m_graph_impl_snapshots.insert({ key, impl_type() });
        }
        for (const auto& key : operators.elevation_snapshot_keys())
        {
            m_graph_impl_snapshots.insert({ key, {} });
        }
    }

    // TODO: probably better to access hydrologically corrected topographic elevation
    // via an explicit flow_graph API method?

    void update_routes()
    {
        impl_type graph_impl;
        impl_type::data_array_type elevation;

        for (const auto& op : m_operators)
        {
            op.apply(graph_impl, elevation);
            op.save(graph_impl, m_graph_impl_snapshots, elevation, m_elevation_snapshots);
        }
    }

    data_array_type m_hydro_elevation;
    graph_impl_map m_graph_impl_snapshots;
    elevation_map m_elevation_snapshots;
    fs::flow_operator_sequence<impl_type> m_operators;
};


void
test()
{
    auto fg = flow_graph_test({ fs::op1(), fs::flow_snapshot("s1"), fs::op2() });
    fg.update_routes();
}


void
add_flow_graph_bindings(py::module& m)
{
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
    py::class_<flow_graph_test>(m, "FlowGraphTest")
        .def(py::init(
            [](const py::list& ops)
            {
                using fg_impl_type = flow_graph_test::impl_type;
                auto op_sequence = fs::make_flow_operator_sequence<fg_impl_type>(ops);
                return std::make_unique<flow_graph_test>(std::move(op_sequence));
            }))
        .def("update_routes", &flow_graph_test::update_routes);

    m.def("test", &test);


    py::class_<fs::py_flow_graph_impl>(m, "FlowGraphImpl")
        .def_property_readonly("receivers", &fs::py_flow_graph_impl::receivers)
        .def_property_readonly("receivers_count", &fs::py_flow_graph_impl::receivers_count)
        .def_property_readonly("receivers_distance", &fs::py_flow_graph_impl::receivers_distance)
        .def_property_readonly("receivers_weight", &fs::py_flow_graph_impl::receivers_weight)
        .def_property_readonly("donors", &fs::py_flow_graph_impl::donors)
        .def_property_readonly("donors_count", &fs::py_flow_graph_impl::donors_count)
        .def_property_readonly("dfs_indices", &fs::py_flow_graph_impl::dfs_indices)
        .def_property_readonly("basins", &fs::py_flow_graph_impl::basins);

    py::class_<fs::py_flow_graph> pyfgraph(m, "FlowGraph");

    fs::add_init_methods<fs::py_profile_grid>(pyfgraph);
    fs::add_init_methods<fs::py_raster_grid>(pyfgraph);

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
