/**
 * @file
 * @brief Provides implementation for the classical Union-Find data structure
 *
 * @author Guillaume Cordonnier
 */

#pragma once

#include "utils.hpp"
#include <vector>

/**
 * @class
 * @brief Union-Find data structure
 *
 * The union find is a standart data structure for storing and merging equivalence
 * classes. Initialy (after construction, or call to clear), the Union Find stores
 * size() distincts classes, denoted by an integer of type T_Class.
 * The classes can be merged with a call to merge(), and obtained through the call
 * to find().
 * Amortized complexity of Union and Find: optimal O(alpha(m, n))
 */
template <class T_Class>
class UnionFind_T
{
public:
    /**
     * @brief UnionFind Constructor
     */
    UnionFind_T()
    {
    }

    /**
     * @brief UnionFind Constructor
     * @param _size the number of classes
     * Complexity O(_size)
     */
    UnionFind_T(size_t _size)
    {
        resize(_size);
    }

    /**
     * @brief clear
     * Restore the union-find data structure. The structure holds m_size different size()
     * distinct classes
     * Complexity O(size())
     */
    void clear()
    {
        size_t old_size = size();
        parent.clear();
        rank.clear();
        resize(old_size);
    }

    /**
     * @brief reserve
     * Allocate some memory for a given number of class
     * @param _size the number of classes
     */
    void reserve(size_t _size)
    {
        parent.reserve(_size);
        rank.reserve(_size);
    }

    /**
     * @brief push_back
     * append a new item at the end of the union find structure
     * @param c the class of the new item
     */

    void push_back(T_Class c)
    {
        parent.push_back(c);
        rank.push_back(0);
        if (c != parent.size() - 1 && !rank[c])
            rank[c] = 1;
    }

    /**
     * @brief resize
     * Resize the internal containers of union find. This only add classes if the new
     * size is larger. Else, some elements are removed, but it is possible that the
     * class returned by find is bigger than the new size.
     * Complexity O(_size - size())
     * @param _size the new size
     */
    void resize(index_t _size)
    {
        size_t old_size = size();

        parent.resize(_size);
        rank.resize(_size, 0);

        std::iota(parent.begin() + old_size, parent.end(), old_size);
    }

    /**
     * @brief size
     * @return the initial number of elements in the union find structure
     */
    size_t size()
    {
        return parent.size();
    }

    /**
     * @brief merge two equivalence class
     * The new equivelence class can be represented either by x or y
     * @param x class to be merged
     * @param y class to be merged
     */
    void merge(T_Class x, T_Class y)
    {
        x = find(x);
        y = find(y);

        if (x != y)
        {
            if (rank[x] < rank[y])
                parent[x] = y;
            else
            {
                parent[y] = x;
                if (rank[x] == rank[y])
                    rank[x] += 1;
            }
        }
    }

    /**
     * @brief find the equivalence class of an element
     * @param x the element for which the class is needed
     * @return the class representative
     */
    T_Class find(T_Class x)
    {
        // find class
        T_Class c = x;
        while (c != parent[c])
            c = parent[c];

        // set class
        while (x != parent[x])
        {
            T_Class t = parent[x];
            parent[x] = c;
            x = t;
        }

        // return class
        return c;
    }

private:
    std::vector<T_Class> parent;
    std::vector<index_t> rank;
};

using UnionFind = UnionFind_T<index_t>;
