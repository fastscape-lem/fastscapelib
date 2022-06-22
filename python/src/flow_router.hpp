#ifndef PYFASTSCAPELIB_FLOW_ROUTER_H
#define PYFASTSCAPELIB_FLOW_ROUTER_H

#include <memory>
#include <stdexcept>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "pytensor_utils.hpp"

#include "fastscapelib/flow/flow_router.hpp"
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
            // corresponding flow_router class)
            virtual ~py_flow_router() = default;

            py_flow_router() = default;
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
            if (dynamic_cast<const py_single_flow_router*>(&py_router))
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
    }
}

#endif
