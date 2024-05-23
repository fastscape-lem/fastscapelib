#ifndef FASTSCAPELIB_FLOW_IMPL_FLOW_GRAPH_INL_HPP
#define FASTSCAPELIB_FLOW_IMPL_FLOW_GRAPH_INL_HPP

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
    template <class G, class S, class Tag>
    flow_graph<G, S, Tag>::flow_graph(G& grid, operators_type operators)
        : m_grid(grid)
        , m_operators(std::move(operators))
        , m_thread_pool(10)
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

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    const std::vector<const flow_operator*>& flow_graph<G, S, Tag>::operators() const
    {
        return m_operators.m_op_vec;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    bool flow_graph<G, S, Tag>::single_flow() const
    {
        return m_operators.out_flowdir() == flow_direction::single;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    const std::vector<std::string>& flow_graph<G, S, Tag>::graph_snapshot_keys() const
    {
        return m_operators.graph_snapshot_keys();
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::graph_snapshot(std::string name) const -> self_type&
    {
        return *(m_graph_snapshots.at(name));
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    const std::vector<std::string>& flow_graph<G, S, Tag>::elevation_snapshot_keys() const
    {
        return m_operators.elevation_snapshot_keys();
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::elevation_snapshot(std::string name) const -> const data_array_type&
    {
        return *(m_elevation_snapshots.at(name));
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::update_routes(const data_array_type& elevation)
        -> const data_array_type&
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
            op->save(*m_impl_ptr, m_graph_impl_snapshots, *elevation_ptr, m_elevation_snapshots);
        }

        return *elevation_ptr;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::grid() const -> grid_type&
    {
        return m_grid;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::size() const -> size_type
    {
        return m_grid.size();
    }

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::grid_shape() const -> shape_type
    {
        // grid shape may have a different type (e.g., from xtensor containers)
        auto shape = m_grid.shape();
        shape_type data_array_shape(shape.begin(), shape.end());
        return data_array_shape;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::impl() const -> const impl_type&
    {
        return *m_impl_ptr;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::impl_ptr() const -> std::shared_ptr<impl_type>
    {
        return m_impl_ptr;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::base_levels() const -> std::vector<size_type>
    {
        const auto& impl_levels = m_impl_ptr->base_levels();
        std::vector<size_type> indices(impl_levels.begin(), impl_levels.end());
        return indices;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    template <class C>
    void flow_graph<G, S, Tag>::set_base_levels(C&& levels)
    {
        if (!m_writeable)
        {
            throw std::runtime_error("cannot set base levels (graph is read-only)");
        }

        m_impl_ptr->set_base_levels(levels);
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::mask() const -> xt_array_t<xt_selector, bool>
    {
        return m_impl_ptr->mask();
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    template <class C>
    void flow_graph<G, S, Tag>::set_mask(C&& mask)
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

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::accumulate(const data_array_type& src) const -> data_array_type
    {
        return m_impl_ptr->accumulate(src);
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    void flow_graph<G, S, Tag>::accumulate(data_array_type& acc, const data_array_type& src) const
    {
        return m_impl_ptr->accumulate(acc, src);
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::accumulate(data_type src) const -> data_array_type
    {
        return m_impl_ptr->accumulate(src);
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    void flow_graph<G, S, Tag>::accumulate(data_array_type& acc, data_type src) const
    {
        return m_impl_ptr->accumulate(acc, src);
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    auto flow_graph<G, S, Tag>::basins() -> data_array_size_type
    {
        data_array_size_type basins = data_array_size_type::from_shape(m_grid.shape());
        auto basins_flat = xt::flatten(basins);

        m_impl_ptr->compute_basins();
        basins_flat = m_impl_ptr->basins();

        return basins;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    int flow_graph<G, S, Tag>::apply_kernel_seq(const FSFlowKernel& kernel, double dt)
    {
        void* new_node_data = kernel.node_data_create();

        for (std::size_t i : impl().dfs_indices())
        {
            if (kernel.node_data_getter(i, kernel.data, new_node_data))
            {
                throw std::runtime_error("Invalid index encountered in node_data getter "
                                         "function\n"
                                         "Please check if you are using dynamic receivers count "
                                         "('max_receivers=-1') or adjust this setting in the "
                                         "'Kernel' "
                                         "specification");
            };
            kernel.func(new_node_data, dt);
            kernel.node_data_setter(i, new_node_data, kernel.data);
        }

        kernel.node_data_free(new_node_data);
        return 0;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    int flow_graph<G, S, Tag>::apply_kernel_seq2(NumbaFlowKernel& kernel, double dt)
    {
        NumbaJitClass new_node_data = kernel.node_data_create();

        for (std::size_t i : impl().dfs_indices())
        {
            if (kernel.node_data_getter(i, kernel.data, new_node_data))
            {
                throw std::runtime_error("Invalid index encountered in node_data getter "
                                         "function\n"
                                         "Please check if you are using dynamic receivers count "
                                         "('max_receivers=-1') or adjust this setting in the "
                                         "'Kernel' "
                                         "specification");
            };
            kernel.func(new_node_data, dt);
            kernel.node_data_setter(i, new_node_data, kernel.data);
        }

        kernel.node_data_free(new_node_data.meminfoptr);
        return 0;
    }

    // int
    // apply_kernel_par(std::vector<std::size_t>& indices, const FSFlowKernel& kernel, double
    // dt)
    // {
    //     std::vector<std::thread> threads;
    //     const blocks<std::size_t> blks(0, indices.size(), kernel.n_threads);

    //     for (std::size_t blk = 0; blk < blks.get_num_blocks(); ++blk)
    //     {
    //         auto run_block = [&kernel, &dt, start = blks.start(blk), end = blks.end(blk)]()
    //         {
    //             void* new_node_data = kernel.node_data_create();

    //             for (std::size_t i = start; i < end; ++i)
    //             {
    //                 if (kernel.node_data_getter(i, kernel.data, new_node_data))
    //                 {
    //                     throw std::runtime_error(
    //                         "Invalid index encountered in node_data getter "
    //                         "function\n"
    //                         "Please check if you are using dynamic receivers count "
    //                         "('max_receivers=-1') or adjust this setting in the "
    //                         "'Kernel' "
    //                         "specification");
    //                 };
    //                 kernel.func(new_node_data, dt);
    //                 kernel.node_data_setter(i, new_node_data, kernel.data);
    //             }

    //             kernel.node_data_free(new_node_data);
    //         };

    //         threads.emplace_back(std::thread(run_block));
    //     }

    //     for (auto& thread : threads)
    //         thread.join();

    //     return 0;
    // }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    int flow_graph<G, S, Tag>::apply_kernel_par2(NumbaFlowKernel& kernel, double dt)
    {
        using bfs_indices_type = typename impl_type::bfs_indices_type;
        const bfs_indices_type *indices, *levels;

        switch (kernel.application_order)
        {
            case kernel_application_order::ANY:
                indices = &impl().storage_indices();
                levels = &impl().random_levels();
                break;
            case kernel_application_order::BREADTH_UPSTREAM:
                indices = &impl().bfs_indices();
                levels = &impl().bfs_levels();
                break;
            default:
                throw std::runtime_error("Unsupported kernel application order");
                break;
        }

        auto n_threads = kernel.n_threads;

        m_thread_pool.resume();
        m_thread_pool.resize(n_threads);

        std::vector<NumbaJitClass> node_data(n_threads);
        for (auto i = 0; i < n_threads; ++i)
            node_data[i] = kernel.node_data_create();

        auto run = [&kernel, &dt, indices, node_data](
                       std::size_t runner, std::size_t start, std::size_t end)
        {
            for (auto i = start; i < end; ++i)
            {
                auto node_idx = (*indices)[i];
                NumbaJitClass n_data = node_data[runner];
                if (kernel.node_data_getter(node_idx, kernel.data, n_data))
                {
                    throw std::runtime_error(
                        "Invalid index encountered in node_data getter "
                        "function\n"
                        "Please check if you are using dynamic receivers count "
                        "('max_receivers=-1') or adjust this setting in the "
                        "'Kernel' "
                        "specification");
                };
                kernel.func(n_data, dt);
                kernel.node_data_setter(node_idx, n_data, kernel.data);
            }
        };

        for (auto i = 1; i < levels->size(); ++i)
        {
            run_blocks((*levels)[i - 1], (*levels)[i], run);
        }

        for (auto i = 0; i < n_threads; ++i)
            kernel.node_data_free(node_data[i].meminfoptr);

        // std::cout << m_thread_pool.was_empty() << std::endl;
        // std::cout << m_thread_pool.started() << std::endl;
        // std::cout << m_thread_pool.stopped() << std::endl;

        m_thread_pool.pause();

        // std::cout << m_thread_pool.paused() << std::endl;

        return 0;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    int flow_graph<G, S, Tag>::apply_kernel(NumbaFlowKernel& kernel, double dt)
    {
        int ret;
        if (kernel.n_threads > 1)
            ret = apply_kernel_par2(kernel, dt);
        else
            ret = apply_kernel_seq2(kernel, dt);

        return ret;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    flow_graph<G, S, Tag>::flow_graph(grid_type& grid, bool single_flow)
        : m_writeable(false)
        , m_grid(grid)
        , m_thread_pool(10)
    {
        m_impl_ptr = std::make_shared<impl_type>(grid, single_flow);
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class G, class S, class Tag>
    template <typename T, typename F, typename R>
    void flow_graph<G, S, Tag>::run_blocks(const T first_index, const T index_after_last, F&& func)
    {
        auto thread_pool_size = m_thread_pool.size();
        std::vector<Job> jobs(thread_pool_size);

        if (index_after_last > first_index)
        {
            const thread_pool::blocks blks(first_index, index_after_last, thread_pool_size);

            for (auto i = 0; i < thread_pool_size; ++i)
            {
                if (i < blks.num_blocks())
                    jobs[i] = Job([i,
                                   func = std::forward<F>(func),
                                   start = blks.start(i),
                                   end = blks.end(i)]() { func(i, start, end); });
                else
                    jobs[i] = Job();
            }
            m_thread_pool.set_tasks(jobs);
            m_thread_pool.run_tasks();

            m_thread_pool.wait();
        };
    }
}  // namespace fastscapelib

#endif  // FASTSCAPELIB_FLOW_IMPL_FLOW_GRAPH_INL_HPP
