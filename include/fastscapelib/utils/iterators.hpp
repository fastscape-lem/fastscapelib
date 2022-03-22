#ifndef FASTSCAPELIB_UTILS_ITERATORS_H
#define FASTSCAPELIB_UTILS_ITERATORS_H

#include "xtl/xiterator_base.hpp"


namespace fastscapelib
{

    template <class G>
    struct index_iterator
        : public xtl::xbidirectional_iterator_base<index_iterator<G>,
                                                   typename G::size_type,
                                                   std::ptrdiff_t,
                                                   typename G::size_type*,
                                                   typename G::size_type>
    {
    public:
        using self_type = index_iterator<G>;
        using base_type = xtl::xbidirectional_iterator_base<self_type,
                                                            typename G::size_type,
                                                            std::ptrdiff_t,
                                                            typename G::size_type*,
                                                            typename G::size_type>;

        using value_type = typename base_type::value_type;
        using reference = typename base_type::reference;
        using pointer = typename base_type::pointer;
        using difference_type = typename base_type::difference_type;

        index_iterator() = default;

        index_iterator(G& grid, value_type position = 0)
            : m_idx(position)
            , m_grid(&grid)
        {
        }

        inline self_type& operator++()
        {
            ++m_idx;
            return *this;
        }

        inline self_type& operator--()
        {
            --m_idx;
            return *this;
        }

        inline reference operator*() const
        {
            return m_idx;
        }

        mutable value_type m_idx = 0;

    private:
        G* m_grid;
    };

    template <class G>
    inline bool operator==(const index_iterator<G>& lhs, const index_iterator<G>& rhs)
    {
        return lhs.m_idx == rhs.m_idx;
    }


    template <class G, class F>
    struct filtered_index_iterator
        : public xtl::xbidirectional_iterator_base<filtered_index_iterator<G, F>,
                                                   typename G::size_type>
    {
    public:
        using self_type = filtered_index_iterator<G, F>;
        using base_type = xtl::xbidirectional_iterator_base<self_type, typename G::size_type>;

        using value_type = typename base_type::value_type;
        using reference = typename base_type::reference;
        using pointer = typename base_type::pointer;
        using difference_type = typename base_type::difference_type;

        filtered_index_iterator(G& grid, F filter_functor, value_type position = 0)
            : m_idx(position)
            , m_grid(grid)
            , m_filter_func(filter_functor)
        {
            if ((position == 0) && !m_filter_func(grid, position))
            {
                ++m_idx;
            }
        }

        inline self_type& operator++()
        {
            do
            {
                ++m_idx;
            } while ((!m_filter_func(m_grid, m_idx)) && (m_idx < m_grid.size()));

            return *this;
        }

        inline self_type& operator--()
        {
            do
            {
                --m_idx;
            } while ((!m_filter_func(m_grid, m_idx)) && (m_idx > 0));

            return *this;
        }

        inline reference operator*() const
        {
            return m_idx;
        }

    private:
        mutable value_type m_idx;

        G& m_grid;

        F m_filter_func;
    };


    template <class G, class F>
    inline bool operator==(const filtered_index_iterator<G, F>& lhs,
                           const filtered_index_iterator<G, F>& rhs)
    {
        return lhs.m_idx == rhs.m_idx;
    }


    namespace detail
    {

        template <class G>
        class node_indices_iterator
        {
        public:
            using iterator_type = index_iterator<G>;

            node_indices_iterator(G& grid)
                : m_grid(grid)
            {
            }

            inline iterator_type begin()
            {
                return m_grid.nodes_indices_begin();
            };

            inline iterator_type end()
            {
                return m_grid.nodes_indices_end();
            };

        private:
            G& m_grid;
        };

    }
}
#endif
