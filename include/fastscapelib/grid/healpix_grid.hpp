#ifndef FASTSCAPELIB_HEALPIX_GRID_H_
#define FASTSCAPELIB_HEALPIX_GRID_H_

#include "healpix_cxx/healpix_base.h"
#include <healpix_cxx/healpix_tables.h>

// conflict between healpix xcomplex macro and xtl xcomplex
#undef xcomplex

#include "fastscapelib/grid/base.hpp"

#include <memory>


namespace fastscapelib
{

    template <class S>
    class healpix_grid;

    /**
     * Healpix grid specialized types
     */
    template <class S>
    struct grid_inner_types<healpix_grid<S>>
    {
        static constexpr bool is_structured = false;
        static constexpr bool is_uniform = false;

        using grid_data_type = double;

        using container_selector = S;
        static constexpr std::size_t container_ndims = 1;

        static constexpr uint8_t n_neighbors_max = 8;
        using neighbors_cache_type = neighbors_no_cache<0>;
    };

    /**
     * @brief 2-dimensional grid on the sphere (HEALPix).
     *
     * Fastscapelib grid adapter for a HEALPix (Hierarchical Equal Area
     * isoLatitude Pixelation of a sphere) grid.
     *
     * @tparam S The container selector for data array members.
     */
    template <class S>
    class healpix_grid : public grid<healpix_grid<S>>
    {
    public:
        using self_type = healpix_grid<S>;
        using base_type = grid<self_type>;

        using grid_data_type = typename base_type::grid_data_type;

        using container_selector = typename base_type::container_selector;
        using size_type = typename base_type::size_type;
        using shape_type = typename base_type::shape_type;

        healpix_grid(std::size_t nside);

    protected:
        using healpix_grid_type = T_Healpix_Base<std::size_t>;
        std::unique_ptr<healpix_grid_type> m_healpix_grid_ptr;

        friend class grid<self_type>;
    };


    /**
     * @name Constructors
     */
    /**
     * Creates a new HEALPix grid
     *
     * @param nside number of divisions along the side of a base-resolution HEALPix pixel.
     */
    template <class S>
    healpix_grid<S>::healpix_grid(std::size_t nside)
    {
        m_healpix_grid_ptr
            = std::make_unique<healpix_grid_type>(nside, Healpix_Ordering_Scheme::NEST);
    }
    //@}
}

#endif  // FASTSCAPELIB_HEALPIX_GRID_H_
