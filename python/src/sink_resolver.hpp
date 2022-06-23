#ifndef PYFASTSCAPELIB_SINK_RESOLVER_H
#define PYFASTSCAPELIB_SINK_RESOLVER_H

#include <memory>
#include <stdexcept>

#include "pytensor_utils.hpp"

#include "fastscapelib/flow/sink_resolver.hpp"
#include "fastscapelib/flow/flow_graph.hpp"


namespace py = pybind11;
namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace detail
    {
        /*
         * Sink resolver facade class for the Python bindings.
         * It does not depend on flow_graph.
         */
        class py_sink_resolver
        {
        public:
            // needed in order for dynamic_cast to work
            // (we'll perform downcasting when creating the
            // corresponding sink_resolver class)
            virtual ~py_sink_resolver() = default;

            py_sink_resolver() = default;
        };

        class py_no_sink_resolver : public py_sink_resolver
        {
        public:
            virtual ~py_no_sink_resolver() = default;

            py_no_sink_resolver() = default;
        };

        /*
         * Create an instance of sink_resolver from a instance of py_sink_resolver and
         * return a unique pointer to that instance.
         */
        template <class FG>
        std::unique_ptr<fs::sink_resolver<FG>> make_sink_resolver(py_sink_resolver& py_resolver)
        {
            if (dynamic_cast<const py_no_sink_resolver*>(&py_resolver))
            {
                return std::make_unique<fs::no_sink_resolver<FG>>();
            }
            else
            {
                throw std::runtime_error("Unsupported sink resolver.");
            }
        }
    }
}

#endif
