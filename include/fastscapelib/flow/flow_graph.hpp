#ifndef FASTSCAPELIB_FLOW_FLOW_GRAPH_H
#define FASTSCAPELIB_FLOW_FLOW_GRAPH_H

#include <map>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <string>
#include <vector>

#include "xtensor/xstrided_view.hpp"

#include "fastscapelib/flow/flow_graph_impl.hpp"
#include "fastscapelib/flow/flow_operator.hpp"
#include "fastscapelib/utils/containers.hpp"


namespace fastscapelib
{

    /**
     * Main class used to compute or follow flow routes on
     * the topographic surface.
     *
     * @tparam G The grid type.
     * @tparam S The xtensor container selector for data array members.
     * @tparam Tag The flow graph implementation tag.
     */
    template <class G,
              class S = typename G::container_selector,
              class Tag = flow_graph_fixed_array_tag>
    class flow_graph
    {
    public:
        using self_type = flow_graph<G, S, Tag>;
        using grid_type = G;
        using container_selector = S;
        using impl_type = detail::flow_graph_impl<G, S, Tag>;
        using operators_type = flow_operator_sequence<impl_type>;

        using size_type = typename grid_type::size_type;

        using data_type = typename grid_type::grid_data_type;
        using data_array_type = dynamic_shape_container_t<container_selector, data_type>;
        using shape_type = typename data_array_type::shape_type;
        using data_array_size_type = dynamic_shape_container_t<container_selector, size_type>;

        using graph_map = std::map<std::string, std::unique_ptr<self_type>>;
        using graph_impl_map = std::map<std::string, std::shared_ptr<impl_type>>;
        using elevation_map = std::map<std::string, std::unique_ptr<data_array_type>>;

        /**
         * Flow graph constructor.
         *
         * @param grid The grid object.
         * @param operators The sequence of flow operators to use for updating the graph.
         */
        flow_graph(G& grid, operators_type operators)
            : m_grid(grid)
            , m_operators(std::move(operators))
        {
            m_impl_ptr = std::make_shared<impl_type>(grid, m_operators.all_single_flow());

            // sanity checks
            if (!m_operators.graph_updated())
            {
                throw std::invalid_argument(
                    "must have at least one operator that updates the flow graph");
            }
            if (m_operators.out_flowdir() == flow_direction::undefined)
            {
                throw std::invalid_argument(
                    "must have at least one operator that defines the output flow direction type");
            }

            // initialize default base levels at fixed value nodes
            m_impl_ptr->set_base_levels(m_grid.nodes_indices(node_status::fixed_value));

            // pre-allocate graph and elevation snapshots
            for (const auto& key : m_operators.graph_snapshot_keys())
            {
                bool single_flow = m_operators.snapshot_single_flow(key);
                auto graph = new self_type(grid, single_flow);

                m_graph_snapshots.insert({ key, std::unique_ptr<self_type>(std::move(graph)) });
                m_graph_impl_snapshots.insert({ key, (*m_graph_snapshots.at(key)).m_impl_ptr });
            }
            for (const auto& key : m_operators.elevation_snapshot_keys())
            {
                auto snapshot = data_array_type::from_shape(grid.shape());
                m_elevation_snapshots.insert(
                    { key, std::make_unique<data_array_type>(std::move(snapshot)) });
            }

            // pre-allocate hydrologically corrected elevation
            if (m_operators.elevation_updated())
            {
                m_elevation_copy = xt::empty<data_type>(grid.shape());
            }
        }

        /**
         * Returns a vector of pointers to the flow operators objects used in
         * this graph.
         *
         * Note: the ``flow_operator`` base class has very limited polymorphism
         * (i.e., only its ``name`` method is virtual). It is still possible to
         * downcast the returned pointers (e.g., using the visitor pattern).
         */
        const std::vector<const flow_operator*>& operators() const
        {
            return m_operators.m_op_vec;
        }

        /**
         * Returns ``true`` if the graph has single flow directions.
         */
        bool single_flow() const
        {
            return m_operators.out_flowdir() == flow_direction::single;
        }

        /*
         * Returns the names of the graph snapshots (if any).
         */
        const std::vector<std::string>& graph_snapshot_keys() const
        {
            return m_operators.graph_snapshot_keys();
        }

        /**
         * Graph snapshot getter.
         *
         * @param name Name of the snapshot.
         */
        self_type& graph_snapshot(std::string name) const
        {
            return *(m_graph_snapshots.at(name));
        }

        /**
         * Returns the names of the elevation snapshots.
         */
        const std::vector<std::string>& elevation_snapshot_keys() const
        {
            return m_operators.elevation_snapshot_keys();
        }

        /**
         * Elevation snapshot getter.
         *
         * @param name Name of the snapshot.
         */
        const data_array_type& elevation_snapshot(std::string name) const
        {
            return *(m_elevation_snapshots.at(name));
        }

        /**
         * Update flow routes from the input topographic surface.
         *
         * This applies in chain the flow operators and takes snapshots (if any).
         *
         * @param elevation The input topographic surface elevation.
         *
         * @return Either the input elevation unchanged or the elevation that
         * has been updated in order to route flow accross closed depressions.
         *
         */
        const data_array_type& update_routes(const data_array_type& elevation)
        {
            if (!m_writeable)
            {
                throw std::runtime_error("cannot update routes (graph is read-only)");
            }

            data_array_type* elevation_ptr;

            if (m_operators.elevation_updated())
            {
                // reset and use hydrologically corrected elevation
                m_elevation_copy = elevation;
                elevation_ptr = &m_elevation_copy;
            }
            else
            {
                // pretty safe to remove the const qualifier (shouldn't be updated)
                elevation_ptr = const_cast<data_array_type*>(&elevation);
            }

            // loop over flow operator implementations
            for (auto op = m_operators.impl_begin(); op != m_operators.impl_end(); ++op)
            {
                op->apply(*m_impl_ptr, *elevation_ptr);
                op->save(
                    *m_impl_ptr, m_graph_impl_snapshots, *elevation_ptr, m_elevation_snapshots);
            }

            return *elevation_ptr;
        }

        /**
         * Returns a reference to the corresponding grid object.
         */
        grid_type& grid() const
        {
            return m_grid;
        }

        /**
         * Returns the total number of nodes in the graph (equals to the size of
         * the grid).
         */
        size_type size() const
        {
            return m_grid.size();
        }

        /**
         * Returns the shape of the grid node arrays.
         */
        shape_type grid_shape() const
        {
            // grid shape may have a different type (e.g., from xtensor containers)
            auto shape = m_grid.shape();
            shape_type data_array_shape(shape.begin(), shape.end());
            return data_array_shape;
        }

        /**
         * Returns a reference to the flow graph implementation.
         */
        const impl_type& impl() const
        {
            return *m_impl_ptr;
        }

        /**
         * Returns a shared pointer to the flow graph implementation.
         *
         * While this might be useful in some cases (e.g., Python bindings), in
         * general accessing the implementation through ``impl()`` (const
         * reference) should be preferred.
         */
        std::shared_ptr<impl_type> impl_ptr() const
        {
            return m_impl_ptr;
        }

        /**
         * Return the indices of the base level nodes.
         */
        std::vector<size_type> base_levels() const
        {
            const auto& impl_levels = m_impl_ptr->base_levels();
            std::vector<size_type> indices(impl_levels.begin(), impl_levels.end());
            return indices;
        }

        /**
         * Clear all existing base level nodes and set new ones.
         *
         * @tparam C Any stl-compatible container type.
         * @param levels The indices of the new base level nodes.
         */
        template <class C>
        void set_base_levels(C&& levels)
        {
            if (!m_writeable)
            {
                throw std::runtime_error("cannot set base levels (graph is read-only)");
            }

            m_impl_ptr->set_base_levels(levels);
        }

        /**
         * Return a mask of where elements with a value ``true`` correspond
         * to grid nodes that are not included in the flow graph.
         */
        dynamic_shape_container_t<container_selector, bool> mask() const
        {
            return m_impl_ptr->mask();
        }

        /**
         * Set a new grid mask.
         *
         * @tparam C Any xtensor-compatible container or expression type.
         * @param mask The new mask where elements with a value ``true``
         * correspond to grid nodes that are not included in the flow graph.
         */
        template <class C>
        void set_mask(C&& mask)
        {
            if (!m_writeable)
            {
                throw std::runtime_error("cannot set mask (graph is read-only)");
            }
            if (!xt::same_shape(mask.shape(), m_grid.shape()))
            {
                throw std::runtime_error("cannot set mask (shape mismatch with grid shape)");
            }

            m_impl_ptr->set_mask(std::forward<C>(mask));
        }

        /**
         * Traverse the flow graph in the top->down direction and accumulate
         * locally produced quantities or fluxes.
         *
         * The local quantitites or fluxes (i.e., source) may for example
         * correspond to precipitation, surface water runoff, sediment flux,
         * etc.
         *
         * The accumulated values represent at each node of the graph the
         * integral of the source over the node upslope contributing area.
         *
         * For example, if the source units are meters (height), the units of
         * the output accumulated values are cubic meters (volume).
         *
         * @param src The source, must be given per area unit and have the same
         * shape than ``grid_shape``.
         *
         * @return The output accumulated values (array of same shape than ``grid_shape``).
         */
        data_array_type accumulate(const data_array_type& src) const
        {
            return m_impl_ptr->accumulate(src);
        }

        /**
         * Traverse the flow graph in the top->down direction and accumulate
         * locally produced quantities or fluxes.
         *
         * This version allows reusing a pre-allocated array to store the output.
         *
         * @param acc The output accumulated values.
         * @param src The source, must be given per area unit.
         *
         * Both `acc` and `src` must of the same shape then ``grid_shape``.
         */
        void accumulate(data_array_type& acc, const data_array_type& src) const
        {
            return m_impl_ptr->accumulate(acc, src);
        }

        /**
         * Traverse the flow graph in the top->down direction and accumulate
         * locally produced quantities or fluxes.
         *
         * @param src The spatially uniform source (scalar).
         *
         * @return The output accumulated values (array of same shape than ``grid_shape``).
         */
        data_array_type accumulate(data_type src) const
        {
            return m_impl_ptr->accumulate(src);
        }

        /**
         * Traverse the flow graph in the top->down direction and accumulate
         * locally produced quantities or fluxes.
         *
         * @param acc The output accumulated values.
         * @param src The spatially uniform source (scalar).
         *
         * Both `acc` and `src` must of the same shape then ``grid_shape``.
         */
        void accumulate(data_array_type& acc, data_type src) const
        {
            return m_impl_ptr->accumulate(acc, src);
        }

        /**
         * Delineate catchments (or basins).
         *
         * A catchment is defined by all adjacent nodes that flow towards the
         * same outlet (or pit) graph node.
         *
         * @return Catchment ids (array of same shape than ``grid_shape``).
         *
         * @note
         * Results may be cached.
         * All masked grid nodes have the same assigned catchment id set by the maximum
         * limit of the integer value range. It is therefore preferable to mask the
         * results prior to, e.g., plotting it.
         */
        data_array_size_type basins()
        {
            data_array_size_type basins = data_array_size_type::from_shape(m_grid.shape());
            auto basins_flat = xt::flatten(basins);

            m_impl_ptr->compute_basins();
            basins_flat = m_impl_ptr->basins();

            return basins;
        }

    private:
        bool m_writeable = true;
        grid_type& m_grid;
        std::shared_ptr<impl_type> m_impl_ptr;
        data_array_type m_elevation_copy;

        graph_map m_graph_snapshots;
        graph_impl_map m_graph_impl_snapshots;
        elevation_map m_elevation_snapshots;

        operators_type m_operators;

        // used internally for creating graph snapshots
        explicit flow_graph(grid_type& grid, bool single_flow)
            : m_writeable(false)
            , m_grid(grid)
        {
            m_impl_ptr = std::make_shared<impl_type>(grid, single_flow);
        }
    };
}

#endif
