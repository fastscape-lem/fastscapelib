/**
 * grids (channel, raster, etc.) hold and manage grid/mesh geometry,
 * topology and boundary conditions.
 */
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xtensor.hpp"

#include "fastscapelib/consts.hpp"
#include "fastscapelib/xtensor_utils.hpp"


namespace fastscapelib
{

//*****************
//* Grid boundaries
//*****************

enum class node_status : std::uint8_t
{
    core = 0,
    fixed_value_boundary = 1,
    fixed_gradient_boundary = 2,
    looped_boundary = 3
};

struct edge_status
{
    node_status left = node_status::core;
    node_status right = node_status::core;

    edge_status(node_status status)
        : left(status), right(status)
    {
    }

    edge_status(node_status status_left, node_status status_right)
        : left(status_left), right(status_right)
    {
    }

    edge_status(std::initializer_list<node_status> edges)
    {
        if (edges.size() != 2)
        {
            throw std::invalid_argument("border status list must have 4 elements: "
                                        "{left, right}");
        }

        auto iter = edges.begin();
        left = *iter++;
        right = *iter++;
    }
};

struct border_status
{
    node_status top = node_status::core;
    node_status right = node_status::core;
    node_status bottom = node_status::core;
    node_status left = node_status::core;

    border_status(node_status status)
        : top(status), right(status), bottom(status), left(status)
    {
    }

    border_status(node_status status_top,
                  node_status status_right,
                  node_status status_bottom,
                  node_status status_left)
        : top(status_top), right(status_right), bottom(status_bottom), left(status_left)
    {
    }

    border_status(std::initializer_list<node_status> borders)
    {
        if (borders.size() != 4)
        {
            throw std::invalid_argument("border status list must have 4 elements: "
                                        "{top, right, bottom, left}");
        }

        auto iter = borders.begin();
        top = *iter++;
        right = *iter++;
        bottom = *iter++;
        left = *iter++;
    }
};


//***************
//* Grid elements
//***************

struct raster_node
{
    size_t row;
    size_t col;
    node_status status;
};

struct node
{
    size_t idx;
    node_status status;
};

struct raster_neighbor
{
    size_t flatten_idx;
    size_t row;
    size_t col;
    node_status status;
    double distance;
};

struct neighbor
{
    size_t idx;
    double distance;
    node_status status;
};


//*******************
//* Profile grid (1D)
//*******************

template <class Tag>
class profile_grid_xt
{
public:

    using xt_container_tag = Tag;
    static const std::size_t xt_container_ndims = 1;
    using status_data_type = std::underlying_type_t<fastscapelib::node_status>;
    using status_xt_type = xt_container_t<xt_container_tag, status_data_type,
                                          xt_container_ndims>;
    using neighbor_list = std::vector<neighbor>;

    static constexpr std::array<std::ptrdiff_t, 3> offsets {0, -1, 1};

    profile_grid_xt(std::size_t size,
                    const double spacing,
                    const edge_status& status_at_edges,
                    const std::vector<node>& status_at_nodes = {});

    neighbor_list neighbors(std::size_t idx) const;

    std::size_t size() const noexcept;
    double spacing() const noexcept;
    const status_xt_type& node_status() const;

private:
    std::size_t m_size;
    double m_spacing;

    status_xt_type m_node_status;
    const edge_status& m_status_at_edges;
    bool has_looped_boundaries = false;

    std::vector<neighbor_list> m_all_neighbors;
    void precompute_neighbors();
};


template <class Tag>
inline profile_grid_xt<Tag>::profile_grid_xt(std::size_t size,
                                             const double spacing,
                                             const edge_status& status_at_edges,
                                             const std::vector<node>& status_at_nodes)
    : m_size(size), m_spacing(spacing), m_status_at_edges(status_at_edges)
{
    m_node_status = xt::zeros<status_data_type>({size});
    m_node_status[0] = status_at_edges.left;
    m_node_status[size-1] = status_at_edges.right;

    for (const node& inode : status_at_nodes)
    {
        m_node_status(inode.idx) = inode.status;
    }

    bool left_looped = status_at_edges.left == node_status::looped_boundary;
    bool right_looped = status_at_edges.right == node_status::looped_boundary;

    if (left_looped ^ right_looped)
    {
        throw std::invalid_argument("inconsistent looped boundary status at edges");
    }
    else if (left_looped && right_looped)
    {
        has_looped_boundaries = true;
    }

    precompute_neighbors();
}

template <class Tag>
void profile_grid_xt<Tag>::precompute_neighbors()
{
    m_all_neighbors.resize(m_size);

    for (std::size_t idx=1; idx<m_size-1; ++idx)
    {
        for (std::size_t k=1; k<3; ++k)
        {
            std::size_t nb_idx = idx + offsets[k];
            neighbor nb = {nb_idx, m_spacing, m_node_status[nb_idx]};
            m_all_neighbors[idx].push_back(nb);
        }
    }

    m_all_neighbors[0].push_back({1, m_spacing, m_node_status[1]});
    m_all_neighbors[m_size-1].push_back({m_size-2, m_spacing, m_node_status[m_size-2]});

    if (has_looped_boundaries)
    {
        m_all_neighbors[0].push_front({m_size-1, m_spacing, m_node_status[m_size-1]});
        m_all_neighbors[m_size-1].push_back({0, m_spacing, m_node_status[0]});
    }
}

template <class Tag>
inline std::size_t profile_grid_xt<Tag>::size() const noexcept
{
    return m_size;
}

template <class Tag>
inline double profile_grid_xt<Tag>::spacing() const noexcept
{
    return m_spacing;
}

template <class Tag>
inline auto profile_grid_xt<Tag>::node_status() const
    -> const profile_grid_xt<Tag>::status_xt_type&
{
    return m_node_status;
}

template <class Tag>
inline auto profile_grid_xt<Tag>::neighbors(std::size_t idx) const
    -> profile_grid_xt<Tag>::neighbor_list
{
    return m_all_neighbors[idx];
}


using profile_grid = profile_grid_xt<xtensor_tag>;


//******************
//* Raster grid (2D)
//******************

enum class raster_connect : std::size_t
{
   king = 4,
   queen = 8
};

} // namespace fastscapelib
