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
    flow_operator_impl_base(std::shared_ptr<OP> ptr)
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

    flow_operator_impl(std::shared_ptr<OP> ptr) = delete;

    int execute(FG& graph_impl) const;
};

enum class flow_direction
{
    undefined,
    single,
    multiple
};

struct flow_operator
{
    inline static const std::string name = "";
    static constexpr bool elevation_updated = false;
    static constexpr bool graph_updated = false;
    static const flow_direction in_flowdir = flow_direction::undefined;
    static const flow_direction out_flowdir = flow_direction::undefined;
};

struct op1 : flow_operator
{
    inline static const std::string name = "op1";
    static constexpr bool elevation_updated = true;
    static constexpr bool graph_updated = false;
    static const flow_direction in_flowdir = flow_direction::undefined;
    static const flow_direction out_flowdir = flow_direction::undefined;
    int param = 1;
};

template <class FG>
class flow_operator_impl<FG, op1, fs::detail::flow_graph_fixed_array_tag>
    : public flow_operator_impl_base<FG, op1>
{
public:
    using base_type = flow_operator_impl_base<FG, op1>;

    flow_operator_impl(std::shared_ptr<op1> ptr)
        : base_type(std::move(ptr)){};

    int execute(FG& graph_impl) const
    {
        return this->p_op->param;
    }
};

struct op2 : flow_operator
{
    inline static const std::string name = "op2";
    static constexpr bool elevation_updated = false;
    static constexpr bool graph_updated = true;
    static const flow_direction in_flowdir = flow_direction::undefined;
    static const flow_direction out_flowdir = flow_direction::single;
    int param = 2;
};

template <class FG>
struct flow_operator_impl<FG, op2, fs::detail::flow_graph_fixed_array_tag>
    : public flow_operator_impl_base<FG, op2>
{
    using base_type = flow_operator_impl_base<FG, op2>;

    flow_operator_impl(std::shared_ptr<op2> ptr)
        : base_type(std::move(ptr)){};

    int execute(FG& graph_impl) const
    {
        return this->p_op->param;
    }
};


/*
 * Flow operator implementation facade.
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
    flow_operator_impl_facade(std::shared_ptr<OP> op)
        : p_op_impl_wrapper(std::make_unique<flow_operator_impl_wrapper<OP>>(std::move(op)))
    {
    }

    flow_operator_impl_facade(flow_operator_impl_facade<FG>& op_impl_facade) = delete;
    explicit flow_operator_impl_facade(flow_operator_impl_facade<FG>&& op_impl_facade)
        : p_op_impl_wrapper(std::move(op_impl_facade.p_op_impl_wrapper))
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
        flow_operator_impl_wrapper(std::shared_ptr<OP> op)
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


/*
 * Immutable sequence of flow operators (e.g., flow routers, sink resolvers)
 * that are applied in chain when updating a flow graph.
 *
 */
template <class FG>
class flow_operator_sequence
{
public:
    using impl_type = FG;
    using operator_impl_type = flow_operator_impl_facade<impl_type>;
    using const_iterator_type = typename std::vector<operator_impl_type>::const_iterator;

    template <class... OPs>
    flow_operator_sequence(OPs&&... ops)
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

    // implement move semantics only (entity of flow_graph)
    flow_operator_sequence(flow_operator_sequence<FG>& op_sequence) = delete;
    flow_operator_sequence(flow_operator_sequence<FG>&& op_sequence)
        : m_op_impl_vec(std::move(op_sequence.m_op_impl_vec))
        , m_elevation_updated(op_sequence.elevation_updated())
        , m_graph_updated(op_sequence.graph_updated())
        , m_out_flowdir(op_sequence.out_flowdir())
        , m_all_single_flow(op_sequence.all_single_flow())
    {
    }

    /*
     * STL-compatible iterators for looping over the operators.
     */
    const_iterator_type begin()
    {
        return m_op_impl_vec.cbegin();
    }

    const_iterator_type end()
    {
        return m_op_impl_vec.cend();
    }

    /*
     * Returns true if at least one operator in the sequence
     * updates topographic elevation.
     */
    bool elevation_updated() const
    {
        return m_elevation_updated;
    }

    /*
     * Returns true if at least one operator in the sequence
     * updates the flow graph.
     */
    bool graph_updated() const
    {
        return m_graph_updated;
    }

    /*
     * Returns the flow direction type (single, multiple or undefined)
     * of the final state of the flow graph, after having applied all operators.
     */
    flow_direction out_flowdir() const
    {
        return m_out_flowdir;
    }

    /*
     * Returns true if all intermediate states of the flow graph
     * have single flow directions.
     */
    bool all_single_flow() const
    {
        return m_all_single_flow;
    }

protected:
    std::vector<operator_impl_type> m_op_impl_vec;

    bool m_elevation_updated = false;
    bool m_graph_updated = false;
    flow_direction m_out_flowdir = flow_direction::undefined;
    bool m_all_single_flow = true;

    template <class OP>
    void add_operator(OP&& op)
    {
        auto op_ptr = std::make_shared<OP>(std::forward<OP>(op));
        add_operator(std::move(op_ptr));
    }

    template <class OP>
    void add_operator(std::shared_ptr<OP> op_ptr)
    {
        if (op_ptr->in_flowdir != flow_direction::undefined && op_ptr->in_flowdir != m_out_flowdir)
        {
            throw std::invalid_argument("flow operator " + op_ptr->name
                                        + " has incompatible input flow directions");
        }
        if (op_ptr->elevation_updated)
        {
            m_elevation_updated = true;
        }
        if (op_ptr->graph_updated)
        {
            m_graph_updated = true;

            if (op_ptr->out_flowdir != flow_direction::undefined)
            {
                m_out_flowdir = op_ptr->out_flowdir;

                if (op_ptr->out_flowdir != flow_direction::single)
                {
                    m_all_single_flow = false;
                }
            }
        }

        m_op_impl_vec.push_back(operator_impl_type(std::move(op_ptr)));
    }

    // Only for bindings (in case the variadic templates constructor cannot be used).
    template <class _FG, class OPs>
    friend flow_operator_sequence<_FG> make_flow_operator_sequence(OPs&& ops);
};


/*
 * Only for Python bindings (need to add operators one at time from a py::list object).
 */
template <class FG, class OPs>
flow_operator_sequence<FG>
make_flow_operator_sequence(OPs&& ops)
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
            throw py::type_error("invalid flow operator");
        }
    }

    return std::move(op_sequence);
}


class flow_graph_test
{
public:
    using impl_type = flow_graph_impl_test<fs::detail::flow_graph_fixed_array_tag>;

    flow_graph_test() = default;
    flow_graph_test(flow_operator_sequence<impl_type> op_sequence)
        : m_op_sequence(std::move(op_sequence))
    {
        if (!m_op_sequence.graph_updated())
        {
            throw std::invalid_argument(
                "must have at least one operator that updates the flow graph");
        }
        if (m_op_sequence.out_flowdir() == flow_direction::undefined)
        {
            throw std::invalid_argument(
                "must have at least one operator that defines the flow direction type");
        }
    }

    std::vector<int> execute()
    {
        std::vector<int> results;
        impl_type arg;

        for (const auto& op_impl : m_op_sequence)
        {
            results.push_back(op_impl.execute(arg));
        }

        return results;
    }

    flow_operator_sequence<impl_type> m_op_sequence;
};


std::vector<int>
test()
{
    auto fg = flow_graph_test({ op1(), op2() });
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
                using fg_impl_type = flow_graph_test::impl_type;
                auto op_sequence = make_flow_operator_sequence<fg_impl_type>(ops);
                return std::make_unique<flow_graph_test>(std::move(op_sequence));
            }))
        .def("execute", &flow_graph_test::execute);

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
