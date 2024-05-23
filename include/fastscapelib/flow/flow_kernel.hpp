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

    struct NumbaJitClass
    {
        void* meminfoptr;  ///< meminfoptr
        void* dataptr;     ///< dataptr
    };

    /////////////////////////////////////////////////////////////////////////////////////////

    struct flow_kernel
    {
        std::function<int(void*)> func;
        std::function<int(std::size_t, void*, void*)> node_data_getter;
        std::function<int(std::size_t, void*, void*)> node_data_setter;
        std::function<void*()> node_data_create;
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

    /////////////////////////////////////////////////////////////////////////////////////////

    // struct flow_kernel
    // {
    //     int (*func)(void*, double);
    //     int (*node_data_getter)(std::size_t, void*, void*);
    //     int (*node_data_setter)(std::size_t, void*, void*);
    //     void* (*node_data_create)();
    //     void (*node_data_free)(void*);
    //     void* data;
    //     int n_threads;
    // };

    /////////////////////////////////////////////////////////////////////////////////////////

    struct NumbaFlowKernel
    {
        int (*func)(NumbaJitClass);                                          ///< flow kernel
        int (*node_data_getter)(std::size_t, NumbaJitClass, NumbaJitClass);  ///< node data getter
        int (*node_data_setter)(std::size_t, NumbaJitClass, NumbaJitClass);  ///< node data setter
        NumbaJitClass (*node_data_create)();         ///< node data allocator
        void (*node_data_free)(void*);               ///< node data destructor
        int n_threads;                               ///< threads count
        kernel_application_order application_order;  ///< traversal order for kernel application

        operator flow_kernel()
        {
            flow_kernel kernel;

            kernel.func = [this](void* node) -> int
            { return func(*reinterpret_cast<NumbaJitClass*>(node)); };

            kernel.node_data_getter = [this](std::size_t index, void* data, void* node_data) -> int
            {
                return node_data_getter(index,
                                        *reinterpret_cast<NumbaJitClass*>(data),
                                        *reinterpret_cast<NumbaJitClass*>(node_data));
            };

            kernel.node_data_setter = [this](std::size_t index, void* node_data, void* data) -> int
            {
                return node_data_setter(index,
                                        *reinterpret_cast<NumbaJitClass*>(node_data),
                                        *reinterpret_cast<NumbaJitClass*>(data));
            };

            kernel.node_data_create = [this]() -> void*
            {
                auto* ptr = new NumbaJitClass;
                *ptr = node_data_create();
                return ptr;
            };

            kernel.node_data_free = [this](void* data) -> void
            { node_data_free(reinterpret_cast<NumbaJitClass*>(data)->meminfoptr); };

            kernel.n_threads = n_threads;
            kernel.application_order = application_order;

            return kernel;
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////////

    struct NumbaFlowKernelData
    {
        NumbaJitClass data;  ///< data

        operator flow_kernel_data()
        {
            flow_kernel_data kernel_data;

            kernel_data.data = &data;

            return kernel_data;
        }
    };
}  // namespace fastscapelib

#endif  // FASTSCAPELIB_FLOW_KERNEL_HPP
