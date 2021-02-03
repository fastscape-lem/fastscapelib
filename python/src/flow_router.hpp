#ifndef PYFASTSCAPELIB_FLOW_ROUTER_H
#define PYFASTSCAPELIB_FLOW_ROUTER_H

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "pytensor_utils.hpp"

#include "fastscapelib/flow_router.hpp"
#include "fastscapelib/flow_router_factory.hpp"
#include "fastscapelib/flow_graph.hpp"

#include <stdexcept>


namespace py = pybind11;
namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace detail
    {

        struct flow_router_method
        {
            virtual ~flow_router_method() = default;

            flow_router_method(fs::flow_router_methods method, std::unique_ptr<flow_router_parameters> params)
                : method(method), parameters(std::move(params))
            {}

            fs::flow_router_methods method;
            std::unique_ptr<flow_router_parameters> parameters;
        };

        struct dummy_flow_router_method : public flow_router_method
        {
            virtual ~dummy_flow_router_method() = default;

            dummy_flow_router_method()
                : flow_router_method(fs::flow_router_methods::dummy,
                                     std::make_unique<no_flow_router_parameters>())
            {};
        };

        struct single_flow_router_method : public flow_router_method
        {
            virtual ~single_flow_router_method() = default;

            single_flow_router_method()
                : flow_router_method(fs::flow_router_methods::single,
                                     std::make_unique<no_flow_router_parameters>())
            {};
        };

        struct multiple_flow_router_method : public flow_router_method
        {    
            virtual ~multiple_flow_router_method() = default;

            multiple_flow_router_method(double p1, double p2)
                : flow_router_method(fs::flow_router_methods::multiple,
                                     std::make_unique<multiple_flow_router_parameters>(p1, p2))
            {}
        };

        template <class G>
        void register_flow_routers()
        {
            using flow_graph_type = fs::flow_graph<G, double, pyarray_selector>;
            using factory = fs::detail::flow_router_factory<flow_graph_type>;

            factory::insert(fs::flow_router_methods::dummy, 
                                    [](const flow_router_parameters& /*params*/) -> typename factory::router_ptr_type 
                                    {
                                        return std::make_unique<fs::dummy_flow_router<flow_graph_type>>(); 
                                    });

            factory::insert(fs::flow_router_methods::single, 
                                    [](const flow_router_parameters& /*params*/) -> typename factory::router_ptr_type 
                                    {
                                        return std::make_unique<fs::single_flow_router<flow_graph_type>>(); 
                                    });

            factory::insert(fs::flow_router_methods::multiple, 
                                    [](const flow_router_parameters& params) -> typename factory::router_ptr_type 
                                    {
                                        auto& derived = static_cast<const multiple_flow_router_parameters&>(params);
                                        return std::make_unique<fs::multiple_flow_router<flow_graph_type>>(derived.p1, derived.p2); 
                                    });
        }
    }
}

#endif
