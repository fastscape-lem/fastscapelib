/**
 * @file
 * @brief Some utilities used internally.
 */
#pragma once

#include <cstdint>


// type used for indexing arrays and array sizes.
using index_t = int64_t;


namespace fastscapelib
{

namespace detail
{


/**
 * @brief Return true if a given (row, col) index is in bounds (false
 *        otherwise).
 */
template<class S>
bool in_bounds(const S& shape, index_t row, index_t col)
{
    return (row >= 0 && row < (index_t) shape[0]
       && col >= 0 && col < (index_t) shape[1]);
}

/**
 * @brief Return a pair of the (row, col) coordinate of a node
 *        in a dem of ncols columns
 */
template<class Node_T>
std::pair<Node_T, Node_T> coords(Node_T node, Node_T ncols)
{
    return {
        node / ncols,
        node % ncols
    };
}

/**
 * @brief Return the node index of a node given its (row, col) coordinate
 *        in a dem of ncols columns
 */
template<class Node_T>
Node_T index(Node_T row, Node_T col, Node_T ncols)
{
    return col + row * ncols;
}

/**
 * @brief proxy for flattened view waiting it to be implmented in xtensor
 * Unsafe in the genral case, use it as a replacement for maintance help ...
 */

template <class Array_T>
class Flattened2D
{
public:

    using value_type = typename Array_T::value_type;
    using size_type  = typename Array_T::size_type;

    Flattened2D(Array_T& ref) : _array_ref {ref} { _ncols = ref.shape()[1];}

    typename Array_T::value_type operator()(size_type index)
    {
        auto c = coords(index, _ncols);
        return _array_ref(std::get<0>(c), std::get<1>(c));
    }

    typename Array_T::value_type at(size_type index)
    {
        auto c = coords(index, _ncols);
        return _array_ref.at(std::get<0>(c), std::get<1>(c));
    }

private:
    Array_T& _array_ref;
    size_type _ncols;
}

template<class Array_T, class Index_T>
typename Array_T::value_type& index_proxy(Array_T& array, Index_T elmt)
{
    auto c = coords(elmt, )
}


}  // namespace detail

}  // namespace fastscapelib
