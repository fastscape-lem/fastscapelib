/**
 * structured grid abstract class.
 */
#ifndef FASTSCAPELIB_GRID_STRUCTURED_GRID_H
#define FASTSCAPELIB_GRID_STRUCTURED_GRID_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>
#include <type_traits>
#include <map>
#include <vector>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xfixed.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xindex_view.hpp"
#include "xtensor/xnoalias.hpp"

#include "xtl/xiterator_base.hpp"

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/utils/utils.hpp"
#include "fastscapelib/utils/xtensor_utils.hpp"


namespace fastscapelib
{

    //***************************
    //* Structured grid interface
    //***************************

    /**
     * @class structured_grid
     * @brief Extends the common grid interface for all structured grid types.
     *
     * This class only defines a basic interface for all structured grid types.
     * It does not embed any data member, this responsibility
     * is delegated to the inheriting classes.
     *
     * @tparam G Derived grid type.
     */
    template <class G>
    class structured_grid : public grid<G>
    {
    public:
        using base_type = grid<G>;
        using inner_types = grid_inner_types<G>;

        using shape_type = typename base_type::shape_type;
        using length_type = typename inner_types::length_type;
        using spacing_type = typename inner_types::spacing_type;

        using spacing_t = std::conditional_t<std::is_arithmetic<spacing_type>::value,
                                             spacing_type,
                                             const spacing_type&>;

        spacing_t spacing() const noexcept;

        length_type length() const noexcept;

        shape_type shape() const noexcept;

    protected:
        using grid<G>::grid;
        ~structured_grid() = default;
    };

    /**
     * @name Grid properties
     */
    //@{
    /**
     * Returns the (uniform) spacing between two adjacent grid nodes.
     *
     * Returns a copy of the value for 1-d grids or a constant reference otherwise.
     */
    template <class G>
    inline auto structured_grid<G>::spacing() const noexcept -> spacing_t
    {
        return this->derived_grid().m_spacing;
    }

    /**
     * Returns the length of the grid for all its dimensions.
     */
    template <class G>
    inline auto structured_grid<G>::length() const noexcept -> length_type
    {
        return this->derived_grid().m_length;
    }

    /**
     * Returns the shape of the grid.
     */
    template <class G>
    inline auto structured_grid<G>::shape() const noexcept -> shape_type
    {
        return this->derived_grid().m_shape;
    }

}

#endif
