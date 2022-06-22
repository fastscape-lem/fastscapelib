#ifndef PYFASTSCAPELIB_FLOW_ROUTER_H
#define PYFASTSCAPELIB_FLOW_ROUTER_H

#include <stdexcept>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "pytensor_utils.hpp"

#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/flow/flow_router_factory.hpp"
#include "fastscapelib/flow/flow_graph.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;


namespace fastscapelib
{
    namespace detail
    {
        /*
         * Flow router facade class for the Python bindings.
         * It does not depend on flow_graph.
         */
        class py_flow_router
        {
        public:
            // needed in order for dynamic_cast to work
            // (we'll perform downcasting when creating the
            // correspondingflow_router class)
            virtual ~py_flow_router() = default;

            py_flow_router() = default;
        };

        class py_dummy_flow_router : public py_flow_router
        {
        public:
            virtual ~py_dummy_flow_router() = default;

            py_dummy_flow_router() = default;
        };

        class py_single_flow_router : public py_flow_router
        {
        public:
            virtual ~py_single_flow_router() = default;

            py_single_flow_router() = default;
        };

        class py_multiple_flow_router : public py_flow_router
        {
        public:
            virtual ~py_multiple_flow_router() = default;

            py_multiple_flow_router(double p1, double p2)
                : m_p1(p1)
                , m_p2(p2)
            {
            }

            double m_p1, m_p2;
        };

        /*
         * Create an instance of flow_router from a instance of py_flow_router and
         * return a unique pointer to that instance.
         */
        template <class FG>
        std::unique_ptr<fs::flow_router<FG>> make_flow_router(py_flow_router& py_router)
        {
            if (dynamic_cast<const py_dummy_flow_router*>(&py_router))
            {
                return std::make_unique<fs::dummy_flow_router<FG>>();
            }
            else if (dynamic_cast<const py_single_flow_router*>(&py_router))
            {
                return std::make_unique<fs::single_flow_router<FG>>();
            }
            else if (const auto* pyr = dynamic_cast<const py_multiple_flow_router*>(&py_router))
            {
                return std::make_unique<fs::multiple_flow_router<FG>>(pyr->m_p1, pyr->m_p2);
            }
            else
            {
                throw std::runtime_error("Unsupported flow router.");
            }
        }


        struct flow_router_method
        {
            virtual ~flow_router_method() = default;

            flow_router_method(fs::flow_router_methods method,
                               std::unique_ptr<flow_router_parameters> params)
                : method(method)
                , parameters(std::move(params))
            {
            }

            fs::flow_router_methods method;
            std::unique_ptr<flow_router_parameters> parameters;
        };

        struct dummy_flow_router_method : public flow_router_method
        {
            virtual ~dummy_flow_router_method() = default;

            dummy_flow_router_method()
                : flow_router_method(fs::flow_router_methods::dummy,
                                     std::make_unique<no_flow_router_parameters>()){};
        };

        struct single_flow_router_method : public flow_router_method
        {
            virtual ~single_flow_router_method() = default;

            single_flow_router_method()
                : flow_router_method(fs::flow_router_methods::single,
                                     std::make_unique<no_flow_router_parameters>()){};
        };

        struct multiple_flow_router_method : public flow_router_method
        {
            virtual ~multiple_flow_router_method() = default;

            multiple_flow_router_method(double p1, double p2)
                : flow_router_method(fs::flow_router_methods::multiple,
                                     std::make_unique<multiple_flow_router_parameters>(p1, p2))
            {
            }
        };

        template <class G>
        void register_flow_routers()
        {
            using flow_graph_type = fs::flow_graph<G, py_selector>;
            using factory = fs::detail::flow_router_factory<flow_graph_type>;

            factory::insert(fs::flow_router_methods::dummy,
                            [](const flow_router_parameters& /*params*/) ->
                            typename factory::router_ptr_type
                            { return std::make_unique<fs::dummy_flow_router<flow_graph_type>>(); });

            factory::insert(
                fs::flow_router_methods::single,
                [](const flow_router_parameters& /*params*/) -> typename factory::router_ptr_type
                { return std::make_unique<fs::single_flow_router<flow_graph_type>>(); });

            factory::insert(
                fs::flow_router_methods::multiple,
                [](const flow_router_parameters& params) -> typename factory::router_ptr_type
                {
                    auto& derived = static_cast<const multiple_flow_router_parameters&>(params);
                    return std::make_unique<fs::multiple_flow_router<flow_graph_type>>(derived.p1,
                                                                                       derived.p2);
                });
        }
    }
}

#endif
