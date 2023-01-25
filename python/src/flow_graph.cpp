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
    using data_array_type = xt::pyarray<double>;
};

/*
 * Flow operator implementation base class.
 *
 * It exposes two methods:
 *
 * - `apply`: may update in-place a flow graph implementation and/or
 *   topographic elevation
 * - `save` : only used by the flow_snapshot operator.
 *
 * By default, these two methods do nothing. A flow operator implementation
 * should re-implement only the `apply` method.
 *
 * The methods are not declared virtual (no polymorphism). Template
 * specialization + type erasure is used instead.
 *
 * A flow operator implementation is decoupled from its flow operator
 * corresponding instance (shared pointer). While flow operators may be
 * instantiated outside of the flow_graph class, flow operator implementations
 * are instantiated inside the flow_graph class as they need the
 * (grid-dependent) flow graph implementation and topographic elevation types.
 *
 * @tparam FG The flow graph implementation type
 * @tparam OP The flow operator type
 *
 */
template <class FG, class OP>
class flow_operator_impl_base
{
public:
    using graph_impl_type = FG;
    using data_array_type = typename graph_impl_type::data_array_type;

    int apply(FG& graph_impl) const
    {
        return 0;
    }

    void save(const FG& graph_impl,
              std::map<std::string, FG>& graph_impl_snapshots,
              const data_array_type& elevation,
              std::map<std::string, data_array_type>& elevation_snapshots) const
    {
    }

protected:
    flow_operator_impl_base(std::shared_ptr<OP> ptr)
        : p_op(std::move(ptr))
    {
    }

    ~flow_operator_impl_base() = default;

    std::shared_ptr<const OP> p_op;
};

/*
 * Flow operator implementation.
 *
 * This template class is not directly constructible. Template specialized
 * implementation classes must be provided for each operator (OP) for at least
 * one of the available flow graph implementations (Tag).
 *
 * @tparam FG The flow graph implementation type (grid-dependent)
 * @tparam OP The flow operator type
 * @tparam Tag The flow graph implementation tag
 *
 */
template <class FG, class OP, class Tag>
class flow_operator_impl : public flow_operator_impl_base<FG, OP>
{
public:
    flow_operator_impl(std::shared_ptr<OP> ptr) = delete;
};

enum class flow_direction
{
    undefined,
    single,
    multiple
};

/*
 * Flow operator.
 *
 * It represents a logical unit that can read and/or modify in-place the flow
 * graph and topographic elevation. It is a common type for flow routers, sink
 * resolvers and snapshots (save intermediate graph and/or elevation states).
 *
 * A flow operator may have one or more implementations, each relative to a
 * specific flow graph implementation. Note: this class and its derived classes
 * only contain the operators' (static) properties and parameters
 * (implementations are in separate classes).
 *
 * Do not use this class directly (it has no implementation). Use instead its
 * derived classes.
 *
 */
class flow_operator
{
public:
    inline static const std::string name = "";
    static constexpr bool elevation_updated = false;
    static constexpr bool graph_updated = false;
    static const flow_direction in_flowdir = flow_direction::undefined;
    static const flow_direction out_flowdir = flow_direction::undefined;
};

class op1 : public flow_operator
{
public:
    inline static const std::string name = "op1";
    static constexpr bool elevation_updated = true;

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

    int apply(FG& graph_impl) const
    {
        return this->p_op->param;
    }
};

class op2 : public flow_operator
{
public:
    inline static const std::string name = "op2";
    static constexpr bool graph_updated = true;
    static const flow_direction out_flowdir = flow_direction::single;

    int param = 2;
};

template <class FG>
class flow_operator_impl<FG, op2, fs::detail::flow_graph_fixed_array_tag>
    : public flow_operator_impl_base<FG, op2>
{
public:
    using base_type = flow_operator_impl_base<FG, op2>;

    flow_operator_impl(std::shared_ptr<op2> ptr)
        : base_type(std::move(ptr)){};

    int apply(FG& graph_impl) const
    {
        return this->p_op->param;
    }
};

/*
 * Flow snapshot operator.
 *
 * A special flow operator used to save intermediate states
 * of the flow graph and/or topographic elevation values while
 * applying the other operators in chain.
 *
 * Those saved states are accessible from the flow_graph object.
 *
 */
class flow_snapshot : public flow_operator
{
public:
    inline static const std::string name = "flow_snapshot";

    flow_snapshot(std::string snapshot_name, bool save_graph = true, bool save_elevation = false)
        : m_snapshot_name(snapshot_name)
        , m_save_graph(save_graph)
        , m_save_elevation(save_elevation)
    {
    }

    const std::string& snapshot_name() const noexcept
    {
        return m_snapshot_name;
    }

    bool save_graph() const noexcept
    {
        return m_save_graph;
    }

    bool save_elevation() const noexcept
    {
        return m_save_elevation;
    }

private:
    std::string m_snapshot_name;
    bool m_save_graph;
    bool m_save_elevation;
};

/*
 * Flow snapshot operator implementation.
 */
template <class FG>
class flow_operator_impl<FG, flow_snapshot, fs::detail::flow_graph_fixed_array_tag>
    : public flow_operator_impl_base<FG, flow_snapshot>
{
public:
    using base_type = flow_operator_impl_base<FG, flow_snapshot>;
    using data_array_type = typename base_type::data_array_type;

    using graph_impl_map = std::map<std::string, FG>;
    using elevation_map = std::map<std::string, data_array_type>;

    flow_operator_impl(std::shared_ptr<flow_snapshot> ptr)
        : base_type(std::move(ptr)){};

    void save(const FG& graph_impl,
              graph_impl_map& graph_impl_snapshots,
              const data_array_type& elevation,
              elevation_map& elevation_snapshots) const
    {
        if (this->p_op->save_graph())
        {
            _save(graph_impl, get_snapshot(graph_impl_snapshots));
        }
        if (this->p_op->save_elevation())
        {
            _save(elevation, get_snapshot(elevation_snapshots));
        }
    }

private:
    FG& get_snapshot(graph_impl_map& graph_impl_snapshots) const
    {
        return graph_impl_snapshots.at(this->p_op->snapshot_name());
    }

    data_array_type& get_snapshot(elevation_map& elevation_snapshots) const
    {
        return elevation_snapshots.at(this->p_op->snapshot_name());
    }

    void _save(const FG& graph_impl, FG& graph_impl_snapshot) const
    {
        py::print("saving graph");
    }

    void _save(const data_array_type& elevation, data_array_type& elevation_snapshot) const
    {
        py::print("saving elevation");
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
    using data_array_type = typename FG::data_array_type;
    using graph_impl_map = std::map<std::string, FG>;
    using elevation_map = std::map<std::string, data_array_type>;

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

    int apply(FG& arg) const
    {
        return p_op_impl_wrapper->apply(arg);
    }

    void save(const FG& graph_impl,
              graph_impl_map& graph_impl_snapshots,
              const data_array_type& elevation,
              elevation_map& elevation_snapshots) const
    {
        p_op_impl_wrapper->save(graph_impl, graph_impl_snapshots, elevation, elevation_snapshots);
    }

    struct flow_operator_impl_wrapper_base
    {
        virtual ~flow_operator_impl_wrapper_base()
        {
        }
        virtual int apply(FG& arg) const = 0;
        virtual void save(const FG& graph_impl,
                          graph_impl_map& graph_impl_snapshots,
                          const data_array_type& elevation,
                          elevation_map& elevation_snapshots) const
            = 0;
    };

    template <class OP>
    class flow_operator_impl_wrapper : public flow_operator_impl_wrapper_base
    {
    public:
        flow_operator_impl_wrapper(std::shared_ptr<OP> op)
            : m_op_impl(std::move(op))
        {
        }

        int apply(FG& arg) const override
        {
            return m_op_impl.apply(arg);
        }

        void save(const FG& graph_impl,
                  graph_impl_map& graph_impl_snapshots,
                  const data_array_type& elevation,
                  elevation_map& elevation_snapshots) const override
        {
            m_op_impl.save(graph_impl, graph_impl_snapshots, elevation, elevation_snapshots);
        }

    private:
        flow_operator_impl<FG, OP, typename FG::impl_tag> m_op_impl;
    };

    std::unique_ptr<flow_operator_impl_wrapper_base> p_op_impl_wrapper;
};


/*
 * Immutable sequence of flow operators (e.g., flow routers, sink resolvers,
 * flow snapshots) that are applied in chain when updating a flow graph.
 *
 * More precisely, it is a container of flow operator implementation (facade)
 * instances.
 *
 * This class is not intended to be used as a stand-alone container. It is used
 * as an entity of flow_graph and can (should) be created implicitly in the
 * flow_graph constructor.
 *
 * @tparam FG The flow graph implementation type
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
        , m_graph_snapshot_keys(op_sequence.graph_snapshot_keys())
        , m_elevation_snapshot_keys(op_sequence.elevation_snapshot_keys())
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

    const std::vector<std::string>& graph_snapshot_keys() const
    {
        return m_graph_snapshot_keys;
    }

    const std::vector<std::string>& elevation_snapshot_keys() const
    {
        return m_elevation_snapshot_keys;
    }

private:
    std::vector<operator_impl_type> m_op_impl_vec;

    std::vector<std::string> m_graph_snapshot_keys;
    std::vector<std::string> m_elevation_snapshot_keys;

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

    // Add an operator of an abitrary type to the sequence and update
    // the sequence properties accordingly.
    //
    // The current (final state) flow direction is updated only if the operator
    // updates the flow graph and explicitly defines an output flow direction
    // type (i.e., single or multiple).
    //
    // Also checks consistency between the current flow direction of the
    // sequence and the expected input flow direction of the operator to add.
    //
    template <class OP>
    void add_operator(std::shared_ptr<OP> op_ptr)
    {
        if constexpr (std::is_same_v<OP, flow_snapshot>)
        {
            update_snapshots(*op_ptr);
        }
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

    void update_snapshots(const flow_snapshot& snapshot)
    {
        const auto& snapshot_name = snapshot.snapshot_name();

        if (snapshot.save_graph())
        {
            m_graph_snapshot_keys.push_back(snapshot_name);
        }
        if (snapshot.save_elevation())
        {
            m_elevation_snapshot_keys.push_back(snapshot_name);
        }
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


class flow_graph_test
{
public:
    using impl_type = flow_graph_impl_test<fs::detail::flow_graph_fixed_array_tag>;
    using data_array_type = impl_type::data_array_type;

    using graph_impl_map = std::map<std::string, impl_type>;
    using elevation_map = std::map<std::string, data_array_type>;

    flow_graph_test() = default;
    flow_graph_test(flow_operator_sequence<impl_type> operators)
        : m_operators(std::move(operators))
    {
        if (!m_operators.graph_updated())
        {
            throw std::invalid_argument(
                "must have at least one operator that updates the flow graph");
        }
        if (m_operators.out_flowdir() == flow_direction::undefined)
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

    std::vector<int> update_routes()
    {
        std::vector<int> results;
        impl_type arg;
        impl_type::data_array_type elevation;

        for (const auto& op : m_operators)
        {
            results.push_back(op.apply(arg));
            op.save(arg, m_graph_impl_snapshots, elevation, m_elevation_snapshots);
        }

        return results;
    }

    graph_impl_map m_graph_impl_snapshots;
    elevation_map m_elevation_snapshots;
    flow_operator_sequence<impl_type> m_operators;
};


std::vector<int>
test()
{
    auto fg = flow_graph_test({ op1(), flow_snapshot("s1"), op2() });
    return fg.update_routes();
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
    py::class_<flow_snapshot, flow_operator, std::shared_ptr<flow_snapshot>>(m, "FlowSnapshot")
        .def(py::init<std::string, bool, bool>(),
             py::arg("snapshot_name"),
             py::arg("save_graph") = true,
             py::arg("save_elevation") = false);
    py::class_<flow_graph_test>(m, "FlowGraphTest")
        .def(py::init(
            [](const py::list& ops)
            {
                using fg_impl_type = flow_graph_test::impl_type;
                auto op_sequence = make_flow_operator_sequence<fg_impl_type>(ops);
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
