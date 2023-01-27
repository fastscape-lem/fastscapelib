#ifndef PYFASTSCAPELIB_FLOW_GRAPH_H
#define PYFASTSCAPELIB_FLOW_GRAPH_H

#include <memory>

#include "pybind11/pybind11.h"

#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_graph_impl.hpp"
#include "fastscapelib/flow/flow_operator.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/flow/flow_snapshot.hpp"
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
            using size_type = std::size_t;

            using donors_type = xt_tensor_t<py_selector, size_type, 2>;
            using donors_count_type = xt_tensor_t<py_selector, neighbors_count_type, 1>;

            using receivers_type = donors_type;
            using receivers_count_type = donors_count_type;
            using receivers_weight_type = xt_tensor_t<py_selector, double, 2>;
            using receivers_distance_type = xt_tensor_t<py_selector, grid_data_type, 2>;

            using dfs_indices_type = xt_tensor_t<py_selector, size_type, 1>;

            using basins_type = xt_tensor_t<py_selector, size_type, 1>;

            virtual ~flow_graph_impl_wrapper_base(){};

            virtual const receivers_type& receivers() const = 0;

            virtual const receivers_count_type& receivers_count() const = 0;

            virtual const receivers_distance_type& receivers_distance() const = 0;

            virtual const receivers_weight_type& receivers_weight() const = 0;

            virtual const donors_type& donors() const = 0;

            virtual const donors_count_type& donors_count() const = 0;

            virtual const dfs_indices_type& dfs_indices() const = 0;

            virtual const basins_type& basins() const = 0;
        };


        template <class FG>
        class flow_graph_impl_wrapper : public flow_graph_impl_wrapper_base
        {
        public:
            using flow_graph_impl_type = FG;

            using size_type = typename flow_graph_impl_wrapper_base::size_type;
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
            using dfs_indices_type = typename flow_graph_impl_wrapper_base::dfs_indices_type;
            using basins_type = typename flow_graph_impl_wrapper_base::basins_type;

            virtual ~flow_graph_impl_wrapper(){};

            flow_graph_impl_wrapper(FG& graph_impl)
                : m_graph_impl(graph_impl){};

            const receivers_type& receivers() const
            {
                return m_graph_impl.receivers();
            };

            const receivers_count_type& receivers_count() const
            {
                return m_graph_impl.receivers_count();
            };

            const receivers_distance_type& receivers_distance() const
            {
                return m_graph_impl.receivers_distance();
            };

            const receivers_weight_type& receivers_weight() const
            {
                return m_graph_impl.receivers_weight();
            };

            const donors_type& donors() const
            {
                return m_graph_impl.donors();
            };

            const donors_count_type& donors_count() const
            {
                return m_graph_impl.donors_count();
            };

            const dfs_indices_type& dfs_indices() const
            {
                return m_graph_impl.dfs_indices();
            };

            const basins_type& basins() const
            {
                return m_graph_impl.basins();
            };

        private:
            flow_graph_impl_type& m_graph_impl;
        };
    }


    class py_flow_graph_impl
    {
    public:
        using size_type = std::size_t;
        using neighbors_count_type = std::uint8_t;
        using grid_data_type = double;

        using donors_type = xt_tensor_t<py_selector, size_type, 2>;
        using donors_count_type = xt_tensor_t<py_selector, neighbors_count_type, 1>;

        using receivers_type = donors_type;
        using receivers_count_type = donors_count_type;
        using receivers_weight_type = xt_tensor_t<py_selector, double, 2>;
        using receivers_distance_type = xt_tensor_t<py_selector, grid_data_type, 2>;

        using dfs_indices_type = xt_tensor_t<py_selector, size_type, 1>;

        using basins_type = xt_tensor_t<py_selector, size_type, 1>;

        template <class FG>
        py_flow_graph_impl(FG& graph_impl)
            : m_wrapper_ptr(std::make_unique<detail::flow_graph_impl_wrapper<FG>>(graph_impl)){};

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

            virtual size_type size() const = 0;
            virtual shape_type grid_shape() const = 0;

            virtual const py_flow_graph_impl& impl() const = 0;

            virtual const data_array_type& update_routes(const data_array_type& elevation) = 0;

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
                m_graph_impl_ptr = std::make_unique<py_flow_graph_impl>(m_graph_ptr->impl());
            }

            virtual ~flow_graph_wrapper(){};

            size_type size() const
            {
                return m_graph_ptr->size();
            };

            shape_type grid_shape() const
            {
                return m_graph_ptr->grid_shape();
            };

            const py_flow_graph_impl& impl() const
            {
                return *m_graph_impl_ptr;
            };

            const data_array_type& update_routes(const data_array_type& elevation)
            {
                return m_graph_ptr->update_routes(elevation);
            };

            void accumulate(data_array_type& acc, const data_array_type& src) const
            {
                return m_graph_ptr->accumulate(acc, src);
            };
            void accumulate(data_array_type& acc, data_type src) const
            {
                return m_graph_ptr->accumulate(acc, src);
            };
            data_array_type accumulate(const data_array_type& src) const
            {
                return m_graph_ptr->accumulate(src);
            };
            data_array_type accumulate(data_type src) const
            {
                return m_graph_ptr->accumulate(src);
            };

            data_array_size_type basins()
            {
                return m_graph_ptr->basins();
            };

        private:
            std::unique_ptr<flow_graph_type> m_graph_ptr;
            std::unique_ptr<py_flow_graph_impl> m_graph_impl_ptr;
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

        template <class G, class O>
        py_flow_graph(G& grid, O operators)
            : m_wrapper_ptr(
                std::make_unique<detail::flow_graph_wrapper<G>>(grid, std::move(operators))){};

        size_type size() const
        {
            return m_wrapper_ptr->size();
        };

        shape_type grid_shape() const
        {
            return m_wrapper_ptr->grid_shape();
        };

        const py_flow_graph_impl& impl() const
        {
            return m_wrapper_ptr->impl();
        };

        const data_array_type& update_routes(const data_array_type& elevation)
        {
            return m_wrapper_ptr->update_routes(elevation);
        };

        void accumulate(data_array_type& acc, const data_array_type& src) const
        {
            return m_wrapper_ptr->accumulate(acc, src);
        };
        void accumulate(data_array_type& acc, data_type src) const
        {
            return m_wrapper_ptr->accumulate(acc, src);
        };
        data_array_type accumulate(const data_array_type& src) const
        {
            return m_wrapper_ptr->accumulate(src);
        };
        data_array_type accumulate(data_type src) const
        {
            return m_wrapper_ptr->accumulate(src);
        };

        data_array_size_type basins()
        {
            return m_wrapper_ptr->basins();
        };

    private:
        std::unique_ptr<detail::flow_graph_wrapper_base> m_wrapper_ptr;
    };


    template <class G>
    void register_py_flow_graph_init(py::class_<py_flow_graph>& pyfg)
    {
        using impl_type =
            typename fs::flow_graph<G, py_selector, flow_graph_fixed_array_tag>::impl_type;

        pyfg.def(py::init(
            [](G& grid, const py::list& ops)
            {
                auto op_sequence = fs::make_flow_operator_sequence<impl_type>(ops);
                return std::make_unique<py_flow_graph>(grid, std::move(op_sequence));
            }));
    }
}

#endif
