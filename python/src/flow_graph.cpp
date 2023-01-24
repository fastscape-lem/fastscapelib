#include <memory>

#include "xtensor-python/pytensor.hpp"

#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/stl.h"

#include "grid.hpp"
#include "flow_graph.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


template <class Tag>
struct flow_graph_impl_test
{
    using impl_tag = Tag;
};


template <class FG, class OP>
class flow_operator_impl_base
{
protected:
    flow_operator_impl_base(std::shared_ptr<OP>&& ptr)
        : p_op(std::move(ptr))
    {
    }

    ~flow_operator_impl_base() = default;

    std::shared_ptr<const OP> p_op;
};


template <class FG, class OP, class Tag>
class flow_operator_impl : public flow_operator_impl_base<FG, OP>
{
public:
    using base_type = flow_operator_impl_base<FG, OP>;

    flow_operator_impl(std::shared_ptr<OP>&& ptr) = delete;

    int execute(FG& graph_impl) const;
};

struct flow_operator
{
};

struct op1 : flow_operator
{
    int param = 1;
};

template <class FG>
class flow_operator_impl<FG, op1, fs::detail::flow_graph_fixed_array_tag>
    : public flow_operator_impl_base<FG, op1>
{
public:
    using base_type = flow_operator_impl_base<FG, op1>;

    flow_operator_impl(std::shared_ptr<op1>&& ptr)
        : base_type(std::move(ptr)){};

    int execute(FG& graph_impl) const
    {
        return this->p_op->param;
    }
};

struct op2 : flow_operator
{
    int param = 2;
};

template <class FG>
struct flow_operator_impl<FG, op2, fs::detail::flow_graph_fixed_array_tag>
    : public flow_operator_impl_base<FG, op2>
{
    using base_type = flow_operator_impl_base<FG, op2>;

    flow_operator_impl(std::shared_ptr<op2>&& ptr)
        : base_type(std::move(ptr)){};

    int execute(FG& graph_impl) const
    {
        return this->p_op->param;
    }
};


/*
 * Flow operator implementation facade, used internally in flow_graph.
 *
 * This facade class implements type erasure so that multiple flow operators
 * of different types can be applied in chain in a flow graph.
 *
 * It takes a flow operator instance as input and creates a wrapper to the
 * corresponding flow operator implementation instance (if such implementation
 * exists for the flow graph implementation given as template argument).
 *
 * @tparam FG The flow graph implementation type
 *
 */
template <class FG>
class flow_operator_impl_facade
{
public:
    template <class OP>
    flow_operator_impl_facade(std::shared_ptr<OP>&& op)
        : p_op_impl_wrapper(std::make_unique<flow_operator_impl_wrapper<OP>>(std::move(op)))
    {
    }

    int execute(FG& arg) const
    {
        return p_op_impl_wrapper->execute(arg);
    }

    struct flow_operator_impl_wrapper_base
    {
        virtual ~flow_operator_impl_wrapper_base()
        {
        }
        virtual int execute(FG& arg) const = 0;
    };

    template <class OP>
    class flow_operator_impl_wrapper : public flow_operator_impl_wrapper_base
    {
    public:
        flow_operator_impl_wrapper(std::shared_ptr<OP>&& op)
            : m_op_impl(std::move(op))
        {
        }

        int execute(FG& arg) const override
        {
            return m_op_impl.execute(arg);
        }

    private:
        flow_operator_impl<FG, OP, typename FG::impl_tag> m_op_impl;
    };

    std::unique_ptr<flow_operator_impl_wrapper_base> p_op_impl_wrapper;
};


class flow_graph_test
{
public:
    using impl_type = flow_graph_impl_test<fs::detail::flow_graph_fixed_array_tag>;
    using operator_impl_type = flow_operator_impl_facade<impl_type>;

    template <class... OPs>
    flow_graph_test(OPs&&... ops)
    {
        int i = 0;
        (
            [&]
            {
                ++i;
                add_operator(std::forward<OPs>(ops));
            }(),
            ...);
    }

    template <class OP>
    void add_operator(OP&& op)
    {
        auto op_ptr = std::make_shared<OP>(std::forward<OP>(op));
        add_operator(std::move(op_ptr));
    }

    template <class OP>
    void add_operator(std::shared_ptr<OP> op_ptr)
    {
        m_op_impl_vec.push_back(operator_impl_type(std::move(op_ptr)));
    }

    std::vector<int> execute()
    {
        std::vector<int> results;
        impl_type arg;

        for (const auto& op_impl : m_op_impl_vec)
        {
            results.push_back(op_impl.execute(arg));
        }

        return results;
    }

    std::vector<operator_impl_type> m_op_impl_vec;
};


std::vector<int>
test()
{
    auto fg = flow_graph_test(op1(), op2());
    return fg.execute();
}


void
add_flow_graph_bindings(py::module& m)
{
    py::class_<flow_operator, std::shared_ptr<flow_operator>>(m, "FlowOperator");
    py::class_<op1, flow_operator, std::shared_ptr<op1>>(m, "Op1")
        .def(py::init<>())
        .def_readwrite("param", &op1::param);
    py::class_<op2, flow_operator, std::shared_ptr<op2>>(m, "Op2")
        .def(py::init<>())
        .def_readwrite("param", &op2::param);
    py::class_<flow_graph_test>(m, "FlowGraphTest")
        .def(py::init(
            [](const py::list& ops)
            {
                auto fgt_ptr = std::make_unique<flow_graph_test>();
                for (auto op : ops)
                {
                    // TODO: smarter way to do this? fold expressions?
                    try
                    {
                        fgt_ptr->add_operator(op.cast<std::shared_ptr<op1>>());
                        continue;
                    }
                    catch (py::cast_error e)
                    {
                    }
                    try
                    {
                        fgt_ptr->add_operator(op.cast<std::shared_ptr<op2>>());
                        continue;
                    }
                    catch (py::cast_error e)
                    {
                        throw py::type_error("invalid flow operator");
                    }
                }

                return fgt_ptr;
            }))
        //.def("add_operator", [](flow_graph_test& f, std::shared_ptr<op1> op_ptr){
        // f.add_operator(op_ptr); }) .def("add_operator", [](flow_graph_test& f,
        // std::shared_ptr<op2> op_ptr){ f.add_operator(op_ptr); })
        .def("execute", &flow_graph_test::execute);

    // m.def("get_param", &get_param);
    // py::class_<flow_ops_sequential>(m, "FlowOpsSequential")
    //     .def(py::init([](std::vector<flow_op>& ops){ return
    //     std::make_unique<flow_ops_sequential>(ops); })) .def_readonly("ops",
    //     &flow_ops_sequential::p_op_wrappers);
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
