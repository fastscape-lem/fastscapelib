#ifndef PYFASTSCAPELIB_FLOW_GRAPH_H
#define PYFASTSCAPELIB_FLOW_GRAPH_H

#include <memory>
#include <variant>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_graph_impl.hpp"
#include "fastscapelib/flow/flow_operator.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/flow/flow_snapshot.hpp"
#include "fastscapelib/flow/sink_resolver.hpp"
#include "fastscapelib/utils/iterators.hpp"
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
            using size_type = std::size_t;
            using grid_data_type = double;

            using donors_type = xt_tensor_t<py_selector, size_type, 2>;
            using donors_count_type = xt_tensor_t<py_selector, size_type, 1>;

            using receivers_type = donors_type;
            using receivers_count_type = donors_count_type;
            using receivers_weight_type = xt_tensor_t<py_selector, double, 2>;
            using receivers_distance_type = xt_tensor_t<py_selector, grid_data_type, 2>;

            using dfs_indices_type = xt_tensor_t<py_selector, size_type, 1>;
            using nodes_indices_iterator_type = stl_container_iterator_wrapper<dfs_indices_type>;

            using basins_type = xt_tensor_t<py_selector, size_type, 1>;

            virtual ~flow_graph_impl_wrapper_base(){};

            virtual bool single_flow() const = 0;

            virtual const receivers_type& receivers() const = 0;

            virtual const receivers_count_type& receivers_count() const = 0;

            virtual const receivers_distance_type& receivers_distance() const = 0;

            virtual const receivers_weight_type& receivers_weight() const = 0;

            virtual const donors_type& donors() const = 0;

            virtual const donors_count_type& donors_count() const = 0;

            virtual const dfs_indices_type& dfs_indices() const = 0;

            virtual nodes_indices_iterator_type nodes_indices_bottomup() const = 0;

            virtual const basins_type& basins() const = 0;
        };


        template <class FG>
        class flow_graph_impl_wrapper : public flow_graph_impl_wrapper_base
        {
        public:
            virtual ~flow_graph_impl_wrapper(){};

            flow_graph_impl_wrapper(const std::shared_ptr<FG>& graph_impl_ptr)
                : m_graph_impl_ptr(graph_impl_ptr)
            {
            }

            bool single_flow() const
            {
                return m_graph_impl_ptr->single_flow();
            }

            const receivers_type& receivers() const
            {
                return m_graph_impl_ptr->receivers();
            };

            const receivers_count_type& receivers_count() const
            {
                return m_graph_impl_ptr->receivers_count();
            };

            const receivers_distance_type& receivers_distance() const
            {
                return m_graph_impl_ptr->receivers_distance();
            };

            const receivers_weight_type& receivers_weight() const
            {
                return m_graph_impl_ptr->receivers_weight();
            };

            const donors_type& donors() const
            {
                return m_graph_impl_ptr->donors();
            };

            const donors_count_type& donors_count() const
            {
                return m_graph_impl_ptr->donors_count();
            };

            const dfs_indices_type& dfs_indices() const
            {
                return m_graph_impl_ptr->dfs_indices();
            };

            nodes_indices_iterator_type nodes_indices_bottomup() const
            {
                return m_graph_impl_ptr->nodes_indices_bottomup();
            }

            const basins_type& basins() const
            {
                return m_graph_impl_ptr->basins();
            };

        private:
            std::shared_ptr<FG> m_graph_impl_ptr;
        };
    }


    class py_flow_graph_impl
    {
    public:
        using size_type = std::size_t;
        using grid_data_type = double;

        using donors_type = xt_tensor_t<py_selector, size_type, 2>;
        using donors_count_type = xt_tensor_t<py_selector, size_type, 1>;

        using receivers_type = donors_type;
        using receivers_count_type = donors_count_type;
        using receivers_weight_type = xt_tensor_t<py_selector, double, 2>;
        using receivers_distance_type = xt_tensor_t<py_selector, grid_data_type, 2>;

        using dfs_indices_type = xt_tensor_t<py_selector, size_type, 1>;
        using nodes_indices_iterator_type = stl_container_iterator_wrapper<dfs_indices_type>;

        using basins_type = xt_tensor_t<py_selector, size_type, 1>;

        template <class FG>
        py_flow_graph_impl(const std::shared_ptr<FG>& graph_impl_ptr)
            : m_wrapper_ptr(
                std::make_unique<detail::flow_graph_impl_wrapper<FG>>(graph_impl_ptr)){};

        bool single_flow() const
        {
            return m_wrapper_ptr->single_flow();
        }

        const receivers_type& receivers() const
        {
            return m_wrapper_ptr->receivers();
        };

        const receivers_count_type& receivers_count() const
        {
            return m_wrapper_ptr->receivers_count();
        };

        const receivers_distance_type& receivers_distance() const
        {
            return m_wrapper_ptr->receivers_distance();
        };

        const receivers_weight_type& receivers_weight() const
        {
            return m_wrapper_ptr->receivers_weight();
        };

        const donors_type& donors() const
        {
            return m_wrapper_ptr->donors();
        };

        const donors_count_type& donors_count() const
        {
            return m_wrapper_ptr->donors_count();
        };

        const dfs_indices_type& dfs_indices() const
        {
            return m_wrapper_ptr->dfs_indices();
        };

        nodes_indices_iterator_type nodes_indices_bottomup() const
        {
            return m_wrapper_ptr->nodes_indices_bottomup();
        }

        const basins_type& basins() const
        {
            return m_wrapper_ptr->basins();
        };

    private:
        std::unique_ptr<detail::flow_graph_impl_wrapper_base> m_wrapper_ptr;
    };


#define TRY_CAST_OPERATOR(TYPE)                                                                    \
    try                                                                                            \
    {                                                                                              \
        op_sequence.add_operator(op.template cast<std::shared_ptr<TYPE>>());                       \
        continue;                                                                                  \
    }                                                                                              \
    catch (py::cast_error e)                                                                       \
    {                                                                                              \
    }

    /*
     * special flow operator sequence constructor for Python bindings
     * (need to add operators one at time from a py::list object).
     */
    template <class FG, class OPs>
    flow_operator_sequence<FG> make_flow_operator_sequence(OPs&& ops)
    {
        flow_operator_sequence<FG> op_sequence;

        for (auto op : ops)
        {
            TRY_CAST_OPERATOR(single_flow_router)
            TRY_CAST_OPERATOR(multi_flow_router)
            TRY_CAST_OPERATOR(pflood_sink_resolver)
            TRY_CAST_OPERATOR(mst_sink_resolver)
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


    /**
     * Flow graph facade class for Python bindings.
     *
     * It implements type erasure in order to expose a single class to Python
     * for all grid types.
     *
     */
    class py_flow_graph;


    namespace detail
    {

        class flow_graph_wrapper_base
        {
        public:
            using size_type = std::size_t;
            using data_type = double;
            using data_array_type = xt_array_t<py_selector, data_type>;
            using shape_type = data_array_type::shape_type;
            using data_array_size_type = xt_array_t<py_selector, size_type>;

            virtual ~flow_graph_wrapper_base(){};

            virtual bool single_flow() const = 0;
            virtual size_type size() const = 0;
            virtual shape_type grid_shape() const = 0;

            virtual const py_flow_graph_impl& impl() const = 0;

            virtual const std::vector<std::string>& graph_snapshot_keys() const = 0;
            virtual std::unique_ptr<flow_graph_wrapper_base> graph_snapshot(std::string name) const
                = 0;
            virtual const std::vector<std::string>& elevation_snapshot_keys() const = 0;
            virtual const data_array_type& elevation_snapshot(std::string name) const = 0;

            virtual const data_array_type& update_routes(const data_array_type& elevation) = 0;

            virtual std::vector<size_type> base_levels() const = 0;
            virtual void set_base_levels(const std::vector<size_type>& levels) = 0;

            virtual xt_array_t<py_selector, bool> mask() const = 0;
            virtual void set_mask(const xt_array_t<py_selector, bool>& mask) = 0;

            virtual void accumulate(data_array_type& acc, const data_array_type& src) const = 0;
            virtual void accumulate(data_array_type& acc, data_type src) const = 0;
            virtual data_array_type accumulate(const data_array_type& src) const = 0;
            virtual data_array_type accumulate(data_type src) const = 0;

            virtual data_array_size_type basins() = 0;
        };

        template <class G>
        class flow_graph_wrapper : public flow_graph_wrapper_base
        {
        public:
            using flow_graph_type = flow_graph<G, py_selector, flow_graph_fixed_array_tag>;
            using operators_type = typename flow_graph_type::operators_type;

            flow_graph_wrapper(G& grid, operators_type operators)
            {
                m_graph_ptr = std::make_unique<flow_graph_type>(grid, std::move(operators));
                m_graph_impl_ptr = std::make_unique<py_flow_graph_impl>(m_graph_ptr->impl_ptr());
            }

            // used for returning graph snapshots
            flow_graph_wrapper(flow_graph_type* snapshot_graph_ptr)
            {
                m_snapshot_graph_ptr = snapshot_graph_ptr;
                m_graph_impl_ptr
                    = std::make_unique<py_flow_graph_impl>(m_snapshot_graph_ptr->impl_ptr());
            }

            virtual ~flow_graph_wrapper(){};

            bool single_flow() const
            {
                return graph().single_flow();
            }

            size_type size() const
            {
                return graph().size();
            }

            shape_type grid_shape() const
            {
                return graph().grid_shape();
            }

            const py_flow_graph_impl& impl() const
            {
                return *m_graph_impl_ptr;
            }

            const std::vector<std::string>& graph_snapshot_keys() const
            {
                return graph().graph_snapshot_keys();
            }

            std::unique_ptr<flow_graph_wrapper_base> graph_snapshot(std::string name) const
            {
                return std::make_unique<flow_graph_wrapper>(&(graph().graph_snapshot(name)));
            }

            const std::vector<std::string>& elevation_snapshot_keys() const
            {
                return graph().elevation_snapshot_keys();
            }

            const data_array_type& elevation_snapshot(std::string name) const
            {
                return graph().elevation_snapshot(name);
            }

            const data_array_type& update_routes(const data_array_type& elevation)
            {
                return graph().update_routes(elevation);
            }

            std::vector<size_type> base_levels() const
            {
                return graph().base_levels();
            }

            void set_base_levels(const std::vector<size_type>& levels)
            {
                graph().set_base_levels(levels);
            }

            xt_array_t<py_selector, bool> mask() const
            {
                return graph().mask();
            }

            void set_mask(const xt_array_t<py_selector, bool>& mask)
            {
                graph().set_mask(mask);
            }

            void accumulate(data_array_type& acc, const data_array_type& src) const
            {
                graph().accumulate(acc, src);
            }
            void accumulate(data_array_type& acc, data_type src) const
            {
                graph().accumulate(acc, src);
            }
            data_array_type accumulate(const data_array_type& src) const
            {
                return graph().accumulate(src);
            }
            data_array_type accumulate(data_type src) const
            {
                return graph().accumulate(src);
            }

            data_array_size_type basins()
            {
                return graph().basins();
            }

        private:
            std::unique_ptr<flow_graph_type> m_graph_ptr;
            std::unique_ptr<py_flow_graph_impl> m_graph_impl_ptr;
            flow_graph_type* m_snapshot_graph_ptr;

            // Allow reusing this wrapper class for both a flow_graph and its snapshot graphs.
            inline flow_graph_type& graph() const
            {
                if (!m_graph_ptr && !m_snapshot_graph_ptr)
                {
                    std::runtime_error("something went wrong (no graph pointer)");
                }

                if (m_graph_ptr)
                {
                    return *m_graph_ptr;
                }
                else
                {
                    return *m_snapshot_graph_ptr;
                }
            }
        };
    }


    class py_flow_graph
    {
    public:
        using size_type = std::size_t;
        using data_type = double;
        using data_array_type = xt_array_t<py_selector, data_type>;
        using shape_type = data_array_type::shape_type;
        using data_array_size_type = xt_array_t<py_selector, size_type>;

        template <class G>
        py_flow_graph(G& grid, const py::list& operators)
            : m_operators(operators)
        {
            using impl_type =
                typename fs::flow_graph<G, py_selector, flow_graph_fixed_array_tag>::impl_type;

            auto op_sequence = fs::make_flow_operator_sequence<impl_type>(operators);
            m_wrapper_ptr
                = std::make_unique<detail::flow_graph_wrapper<G>>(grid, std::move(op_sequence));
        }

        // used for returning graph snapshots
        py_flow_graph(std::unique_ptr<detail::flow_graph_wrapper_base> graph_ptr)
            : m_wrapper_ptr(std::move(graph_ptr))
        {
        }

        py::list operators() const
        {
            return m_operators;
        }

        bool single_flow() const
        {
            return m_wrapper_ptr->single_flow();
        }

        size_type size() const
        {
            return m_wrapper_ptr->size();
        }

        shape_type grid_shape() const
        {
            return m_wrapper_ptr->grid_shape();
        }

        const py_flow_graph_impl& impl() const
        {
            return m_wrapper_ptr->impl();
        }

        const std::vector<std::string>& graph_snapshot_keys() const
        {
            return m_wrapper_ptr->graph_snapshot_keys();
        }

        std::unique_ptr<py_flow_graph> graph_snapshot(std::string name) const
        {
            return std::make_unique<py_flow_graph>(m_wrapper_ptr->graph_snapshot(name));
        }

        const std::vector<std::string>& elevation_snapshot_keys() const
        {
            return m_wrapper_ptr->elevation_snapshot_keys();
        }

        const data_array_type& elevation_snapshot(std::string name) const
        {
            return m_wrapper_ptr->elevation_snapshot(name);
        }

        const data_array_type& update_routes(const data_array_type& elevation)
        {
            return m_wrapper_ptr->update_routes(elevation);
        }

        std::vector<size_type> base_levels()
        {
            return m_wrapper_ptr->base_levels();
        }

        void set_base_levels(const std::vector<size_type>& levels)
        {
            m_wrapper_ptr->set_base_levels(levels);
        }

        xt_array_t<py_selector, bool> mask() const
        {
            return m_wrapper_ptr->mask();
        }

        void set_mask(const xt_array_t<py_selector, bool>& mask)
        {
            m_wrapper_ptr->set_mask(mask);
        }

        void accumulate(data_array_type& acc, const data_array_type& src) const
        {
            m_wrapper_ptr->accumulate(acc, src);
        }
        void accumulate(data_array_type& acc, data_type src) const
        {
            m_wrapper_ptr->accumulate(acc, src);
        }
        data_array_type accumulate(const data_array_type& src) const
        {
            return m_wrapper_ptr->accumulate(src);
        }
        data_array_type accumulate(data_type src) const
        {
            return m_wrapper_ptr->accumulate(src);
        }

        data_array_size_type basins()
        {
            return m_wrapper_ptr->basins();
        }

    private:
        std::unique_ptr<detail::flow_graph_wrapper_base> m_wrapper_ptr;
        py::list m_operators;
    };


    template <class G>
    void register_py_flow_graph_init(py::class_<py_flow_graph>& pyfg, bool add_docstrings = false)
    {
        auto init_func = [](G& grid, const py::list& operators)
        { return std::make_unique<py_flow_graph>(grid, operators); };

        if (add_docstrings)
        {
            pyfg.def(py::init(init_func),
                     R"doc(__init__(self, grid: Any, operators: List[FlowOperator]) -> None

                     FlowGraph initializer.

                     Parameters
                     ----------
                     grid : object
                         Any Fastscapelib grid object.
                     operators : list
                         A list of :class:`FlowOperator` instances that will be
                         applied in chain (in the given order) when computing
                         the graph (flow routing).

                     )doc");
        }
        else
        {
            pyfg.def(py::init(init_func));
        }
    }

    template <class OP>
    void register_operator_static_attrs(
        py::class_<OP, fs::flow_operator, std::shared_ptr<OP>>& pyop)
    {
        pyop.def_readonly_static("elevation_updated",
                                 &OP::elevation_updated,
                                 "True if the operator updates the input elevation")
            .def_readonly_static(
                "graph_updated", &OP::graph_updated, "True if the operator updates the graph")
            .def_readonly_static(
                "in_flowdir", &OP::in_flowdir, "Expected :class:`FlowDirection` of the input graph")
            .def_readonly_static(
                "out_flowdir", &OP::out_flowdir, ":class:`FlowDirection` of the output graph");
    }
}

#endif
