#ifndef PYFASTSCAPELIB_SINK_RESOLVER_H
#define PYFASTSCAPELIB_SINK_RESOLVER_H

#include <memory>
#include <stdexcept>

#include "pytensor_utils.hpp"

#include "fastscapelib/flow/sink_resolver.hpp"
#include "fastscapelib/flow/sink_resolver_factory.hpp"
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

        struct sink_resolver_method
        {
            virtual ~sink_resolver_method() = default;

            sink_resolver_method(fs::sink_resolver_methods method)
                : method(method)
            {
            }

            fs::sink_resolver_methods method;
        };

        struct no_sink_resolver_method : public sink_resolver_method
        {
            virtual ~no_sink_resolver_method() = default;

            no_sink_resolver_method()
                : sink_resolver_method(fs::sink_resolver_methods::none){};
        };

        template <class G>
        void register_sink_resolvers()
        {
            using flow_graph_type = fs::flow_graph<G, py_selector>;
            using factory = fs::detail::sink_resolver_factory<flow_graph_type>;

            factory::insert(fs::sink_resolver_methods::none,
                            []() -> typename factory::resolver_ptr_type
                            { return std::make_unique<fs::no_sink_resolver<flow_graph_type>>(); });
        }
    }
}

#endif
