#ifndef FASTSCAPELIB_FLOW_KERNEL_HPP
#define FASTSCAPELIB_FLOW_KERNEL_HPP

#include <cstddef>
#include <cstdint>
#include <thread>
#include <iostream>
#include <vector>
#include <functional>

namespace fastscapelib
{
    /**
     * Flow graph traversal direction and order
     *
     * This enum class is used to specify the direction and order in which to
     * visit each node of a flow graph.
     *
     */
    enum class flow_graph_traversal_dir
    {
        any,                /**< unspecified direction */
        depth_downstream,   /**< from up to downstream in the depth-first order */
        depth_upstream,     /**< from down to upstream in the depth-first order */
        breadth_downstream, /**< from up to downstream in the breath-first order */
        breadth_upstream    /**< from down to upstream in the breath-first order */
    };

    /////////////////////////////////////////////////////////////////////////////////////////

    struct flow_kernel
    {
        std::function<int(void*)> func;
        std::function<int(std::size_t, void*, void*)> node_data_getter;
        std::function<int(std::size_t, void*, void*)> node_data_setter;
        std::function<void*()> node_data_create;
        std::function<void(void*, void*)> node_data_init;
        std::function<void(void*)> node_data_free;
        int n_threads;
        int min_block_size;
        int min_level_size;
        flow_graph_traversal_dir apply_dir
            = flow_graph_traversal_dir::any;  ///< order for kernel application
    };

    /////////////////////////////////////////////////////////////////////////////////////////

    struct flow_kernel_data
    {
        void* data;
    };
}  // namespace fastscapelib

#endif  // FASTSCAPELIB_FLOW_KERNEL_HPP
