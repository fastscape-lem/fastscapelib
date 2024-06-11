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
    enum class kernel_application_order
    {
        ANY,
        DEPTH_DOWNSTREAM,
        DEPTH_UPSTREAM,
        BREADTH_DOWNSTREAM,
        BREADTH_UPSTREAM
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
        kernel_application_order application_order
            = kernel_application_order::ANY;  ///< order for kernel application
    };

    /////////////////////////////////////////////////////////////////////////////////////////

    struct flow_kernel_data
    {
        void* data;
    };
}  // namespace fastscapelib

#endif  // FASTSCAPELIB_FLOW_KERNEL_HPP
