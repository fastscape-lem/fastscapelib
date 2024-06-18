#ifndef FASTSCAPELIB_HEALPIX_GRID_H_
#define FASTSCAPELIB_HEALPIX_GRID_H_

#include "healpix_cxx/healpix_base.h"
#include "healpix_cxx/healpix_tables.h"

// conflict between healpix xcomplex macro and xtl xcomplex
#undef xcomplex
#include <xtensor/xbroadcast.hpp>

#include "fastscapelib/grid/base.hpp"

#include <math.h>
#include <memory>


namespace fastscapelib
{

    template <class S, class T>
    class healpix_grid;

    /**
     * Healpix grid specialized types
     */
    template <class S, class T>
    struct grid_inner_types<healpix_grid<S, T>>
    {
        static constexpr bool is_structured = false;
        static constexpr bool is_uniform = false;

        using grid_data_type = double;

        using container_selector = S;
        static constexpr std::size_t container_ndims = 1;

        static constexpr uint8_t n_neighbors_max = 8;
        using neighbors_cache_type = neighbors_no_cache<n_neighbors_max>;
    };

    /**
     * @brief 2-dimensional grid on the sphere (HEALPix).
     *
     * Fastscapelib grid adapter for a HEALPix (Hierarchical Equal Area
     * isoLatitude Pixelation of a sphere) grid.
     *
     * @tparam S The container selector for data array members.
     * @tparam T The integer type used to store the HEALPix grid node indices.
     */
    template <class S, class T = int>
    class healpix_grid : public grid<healpix_grid<S, T>>
    {
    public:
        using self_type = healpix_grid<S, T>;
        using base_type = grid<self_type>;
        using inner_types = grid_inner_types<self_type>;

        using grid_data_type = typename base_type::grid_data_type;

        using container_selector = typename base_type::container_selector;
        using container_type = fixed_shape_container_t<container_selector,
                                                       grid_data_type,
                                                       inner_types::container_ndims>;

        using size_type = typename base_type::size_type;
        using shape_type = typename base_type::shape_type;

        healpix_grid(T nside, double radius = 6.371e6);

        T nside() const;

    protected:
        using healpix_type = T_Healpix_Base<T>;
        std::unique_ptr<healpix_type> m_healpix_obj_ptr;

        shape_type m_shape;
        size_type m_size;
        double m_radius;
        double m_node_area;

        inline container_type nodes_areas_impl() const;
        inline grid_data_type nodes_areas_impl(const size_type& idx) const noexcept;

        static constexpr std::size_t dimension_impl() noexcept;

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
    template <class S, class T>
    healpix_grid<S, T>::healpix_grid(T nside, double radius)
        : m_radius(radius)
    {
        m_healpix_obj_ptr = std::make_unique<healpix_type>(nside, Healpix_Ordering_Scheme::NEST);

        m_size = m_healpix_obj_ptr->Npix();
        m_shape = { static_cast<typename shape_type::value_type>(m_size) };
        m_node_area = 4.0 * M_PI * m_radius * m_radius / m_size;
    }
    //@}

    template <class S, class T>
    auto healpix_grid<S, T>::nside() const -> T
    {
        return m_healpix_obj_ptr->NSide();
    }

    template <class S, class T>
    inline auto healpix_grid<S, T>::nodes_areas_impl() const -> container_type
    {
        return xt::broadcast(m_node_area, m_shape);
    }

    template <class S, class T>
    inline auto healpix_grid<S, T>::nodes_areas_impl(const size_type& /*idx*/) const noexcept
        -> grid_data_type
    {
        return m_node_area;
    }

    template <class S, class T>
    constexpr std::size_t healpix_grid<S, T>::dimension_impl() noexcept
    {
        return 2;
    }
}

#endif  // FASTSCAPELIB_HEALPIX_GRID_H_
