#ifndef FASTSCAPELIB_FLOW_SINK_RESOLVER_FACTORY_H
#define FASTSCAPELIB_FLOW_SINK_RESOLVER_FACTORY_H

#include <functional>
#include <map>

#include "fastscapelib/flow/sink_resolver.hpp"


namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace detail
    {

        /**
         * A ``sink_resolver`` factory to register builders
         * and use them to create new ``sink_resolver``s.
         *
         * @tparam FG The flow_graph class.
         */
        template <class FG>
        class sink_resolver_factory
        {
        public:
            using self_type = sink_resolver_factory<FG>;
            using resolver_ptr_type = std::unique_ptr<fs::sink_resolver<FG>>;
            using func_type = std::function<resolver_ptr_type()>;

            sink_resolver_factory(const sink_resolver_factory&) = delete;
            sink_resolver_factory(sink_resolver_factory&&) = delete;
            sink_resolver_factory& operator=(const sink_resolver_factory&) = delete;
            sink_resolver_factory& operator=(sink_resolver_factory&&) = delete;

            static resolver_ptr_type build(fs::sink_resolver_methods method)
            {
                return get_instance().build_impl(method);
            }

            static bool insert(const fs::sink_resolver_methods& method, func_type&& builder)
            {
                return get_instance().insert_impl(method, std::move(builder));
            }

            static self_type& get_instance()
            {
                static self_type instance;
                return instance;
            }

        private:
            sink_resolver_factory() = default;
            ~sink_resolver_factory() = default;

            resolver_ptr_type build_impl(fs::sink_resolver_methods method) const
            {
                auto iter = m_factory.find(method);
                if (iter != m_factory.end())
                {
                    return (iter->second)();
                }
                else
                {
                    return nullptr;
                }
            }

            bool insert_impl(const fs::sink_resolver_methods& method, func_type&& builder)
            {
                if (m_factory.find(method) != m_factory.end())
                {
                    return false;
                }
                m_factory.insert(std::make_pair(method, builder));
                return true;
            }

            std::map<fs::sink_resolver_methods, func_type> m_factory;
        };

    }
}

#endif
