#ifndef PYFASTSCAPELIB_SINK_RESOLVER_H
#define PYFASTSCAPELIB_SINK_RESOLVER_H

#include "pytensor_utils.hpp"

#include "fastscapelib/flow/sink_resolver.hpp"
#include "fastscapelib/flow/sink_resolver_factory.hpp"
#include "fastscapelib/flow/flow_graph.hpp"

#include <stdexcept>


namespace py = pybind11;
namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace detail
    {

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
            using flow_graph_type = fs::flow_graph<G, double, pyarray_selector>;
            using factory = fs::detail::sink_resolver_factory<flow_graph_type>;

            factory::insert(fs::sink_resolver_methods::none,
                            []() -> typename factory::resolver_ptr_type
                            { return std::make_unique<fs::no_sink_resolver<flow_graph_type>>(); });
        }
    }
}

#endif
