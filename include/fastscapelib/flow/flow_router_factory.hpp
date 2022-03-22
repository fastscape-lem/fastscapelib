#ifndef FASTSCAPELIB_FLOW_FLOW_ROUTER_FACTORY_H
#define FASTSCAPELIB_FLOW_FLOW_ROUTER_FACTORY_H

#include "fastscapelib/flow/flow_router.hpp"

#include <functional>


namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace detail
    {

        /**
         * Base class to define a ``flow_router`` parameters.
         */
        struct flow_router_parameters
        {
        };


        /**
         * Class to define a ``flow_router`` parameters.
         */
        struct no_flow_router_parameters : public flow_router_parameters
        {
        };


        /**
         * Class to define a ``multiple_flow_router`` parameters.
         */
        struct multiple_flow_router_parameters : public flow_router_parameters
        {
            multiple_flow_router_parameters(double param1, double param2)
                : p1(param1)
                , p2(param2)
            {
            }

            double p1;
            double p2;
        };


        /**
         * A ``flow_router`` factory to register builders
         * and use them to create new ``flow_router``s.
         *
         * @tparam FG The flow_graph class.
         */
        template <class FG>
        class flow_router_factory
        {
        public:
            using self_type = flow_router_factory<FG>;
            using router_ptr_type = std::unique_ptr<fs::flow_router<FG>>;
            using func_type = std::function<router_ptr_type(const flow_router_parameters&)>;

            flow_router_factory(const flow_router_factory&) = delete;
            flow_router_factory(flow_router_factory&&) = delete;
            flow_router_factory& operator=(const flow_router_factory&) = delete;
            flow_router_factory& operator=(flow_router_factory&&) = delete;

            static router_ptr_type build(fs::flow_router_methods method,
                                         const flow_router_parameters& params)
            {
                return get_instance().build_impl(method, params);
            }

            static bool insert(const fs::flow_router_methods& method, func_type&& builder)
            {
                return get_instance().insert_impl(method, std::move(builder));
            }

            static self_type& get_instance()
            {
                static self_type instance;
                return instance;
            }

        private:
            flow_router_factory() = default;
            ~flow_router_factory() = default;

            router_ptr_type build_impl(fs::flow_router_methods method,
                                       const flow_router_parameters& params) const
            {
                auto iter = m_factory.find(method);
                if (iter != m_factory.end())
                {
                    return (iter->second)(params);
                }
                else
                {
                    return nullptr;
                }
            }

            bool insert_impl(const fs::flow_router_methods& method, func_type&& builder)
            {
                if (m_factory.find(method) != m_factory.end())
                {
                    return false;
                }
                m_factory.insert(std::make_pair(method, builder));
                return true;
            }

            std::map<fs::flow_router_methods, func_type> m_factory;
        };
    }
}

#endif
