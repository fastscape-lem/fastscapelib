#ifndef FASTSCAPELIB_FLOW_BASIN_GRAPH_H
#define FASTSCAPELIB_FLOW_BASIN_GRAPH_H

#pragma once

#include <vector>
#include <array>
#include <queue>
#include <stack>
#include <limits>
#include <numeric>
#include <algorithm>

#include "fastscapelib/grid/base.hpp"
#include "fastscapelib/flow/flow_graph_impl.hpp"
#include "fastscapelib/utils/utils.hpp"
#include "fastscapelib/utils/consts.hpp"
#include "fastscapelib/utils/union_find.hpp"


class BasinGraph_Test;

namespace fastscapelib
{

    enum class mst_method
    {
        kruskal,
        boruvka
    };

    enum class sink_route_method
    {
        basic,
        carve,
        fill_sloped
    };


    template <class FG>
    class basin_graph
    {
    public:
        using flow_graph_impl_type = FG;

        static_assert(std::is_same<typename flow_graph_impl_type::tag,
                                   detail::flow_graph_fixed_array_tag>::value,
                      "basin graph requires the fixed array flow graph implementation");

        using size_type = typename flow_graph_impl_type::size_type;

        using data_type = typename flow_graph_impl_type::data_type;
        using data_array_type = typename flow_graph_impl_type::data_array_type;

        /*
         * Represents an edge of the graph of flow basins. It contains the indices
         * of the two connected basins as well as the indices of the grid nodes
         * forming the pass crossing the basins (+ the elevation of the pass).
         */
        struct edge
        {
            size_type link[2];
            size_type pass[2];
            data_type pass_elevation;

            static edge make_edge(const size_type& from, const size_type& to)
            {
                return edge{ { from, to },
                             { size_type(-1), size_type(-1) },
                             std::numeric_limits<data_type>::lowest() };
            }

            bool operator==(const edge& other)
            {
                return link[0] == other.link[0] && link[1] == other.link[1]
                       && pass[0] == other.pass[0] && pass[1] == other.pass[1]
                       && pass_elevation == other.weight;
            }
        };


        basin_graph(flow_graph_impl_type& flow_graph_impl,
                    mst_method basin_method,
                    sink_route_method route_method)
            : m_flow_graph_impl(flow_graph_impl)
            , m_mst_method(basin_method)
            , m_sink_route_method(route_method)
        {
            m_perf_boruvka = -1;

            // TODO: check shape of receivers (should be single flow)
        }

        size_type basins_count() const
        {
            return m_outlets.size();
        }

        std::vector<size_type>& outlets() const
        {
            return m_outlets;
        }

        const std::vector<edge>& edges() const
        {
            return m_edges;
        }

        const std::vector<size_type>& tree() const
        {
            return m_tree;
        }

        void compute_basins();

        void update_routes(const data_array_type& elevation);

        template <mst_method algo,
                  sink_route_method connect,
                  class Basins_XT,
                  class Rcv_XT,
                  class DistRcv_XT,
                  class Stack_XT,
                  class Active_XT,
                  class Elevation_XT>
        void update_receivers(Rcv_XT& receivers,
                              DistRcv_XT& dist2receivers,
                              const Basins_XT& basins,
                              const Stack_XT& stack,
                              const Active_XT& active_nodes,
                              const Elevation_XT& elevation,
                              data_type dx,
                              data_type dy);

        size_t perf_boruvka() const
        {
            return m_perf_boruvka;
        }

    protected:
        void connect_basins(const data_array_type& elevation);

        void compute_tree_kruskal();

        void compute_tree_boruvka();

        void orient_edges();

        template <class Rcv_XT, class DistRcv_XT, class Elevation_XT>
        void update_pits_receivers(Rcv_XT& receivers,
                                   DistRcv_XT& dist2receivers,
                                   const Elevation_XT& elevation,
                                   double dx,
                                   double dy);

        template <class Rcv_XT, class DistRcv_XT, class Elevation_XT>
        void update_pits_receivers_carve(Rcv_XT& receivers,
                                         DistRcv_XT& dist2receivers,
                                         const Elevation_XT& elevation,
                                         double dx,
                                         double dy);

        template <class Rcv_XT, class DistRcv_XT, class Elevation_XT, class Basins_XT>
        void update_pits_receivers_sloped(Rcv_XT& receivers,
                                          DistRcv_XT& dist2receivers,
                                          const Elevation_XT& elevation,
                                          const Basins_XT&,
                                          double dx,
                                          double dy);

    private:
        flow_graph_impl_type& m_flow_graph_impl;

        mst_method m_mst_method;
        sink_route_method m_sink_route_method;

        std::vector<size_type> m_outlets;  // bottom nodes of basins
        std::vector<edge> m_edges;
        std::vector<size_type> m_tree;  // indices of edges

        // root is used as a virtual basin graph node to which all outer basins are
        // connected, thus collecting the flow from the whole modeled domain. It
        // is used for the computation of the basin tree (minimum spanning
        // tree). As an implementation detail, this virtual node is actually
        // assigned to one of the outer basins, chosen arbitrarily.
        size_type m_root;

        // optimization for basin connections (arrays storing the positions of
        // already created edges)
        std::vector<size_type> m_edge_positions;
        std::vector<size_type> m_edge_positions_tmp;

        // kruskal
        std::vector<size_type> m_edges_indices;
        detail::union_find<size_type> m_basins_uf;

        // boruvka
        std::vector<std::array<size_type, 2>> m_link_basins;

        struct Connect
        {
            size_type begin;
            size_type size;
        };

        std::vector<Connect> m_adjacency;

        struct EdgeParse
        {
            size_type link_id;
            size_type next;
        };

        std::vector<EdgeParse> m_adjacency_list;
        std::vector<size_type> m_low_degrees;
        std::vector<size_type> m_large_degrees;
        std::vector<size_type> m_edge_bucket;
        std::vector<size_type> m_edge_in_bucket;

        // TODO: make it an option
        // 16 for 8-connectivity raster grid, 8 for plannar graph
        int m_max_low_degree = 16;

        template <class T>
        inline void check_capacity(std::vector<T> vec) const
        {
            assert(vec.size() < vec.capacity());
        }

        size_t m_perf_boruvka;

        inline void increase_perf_boruvka()
        {
            ++m_perf_boruvka;
        }

        // reorder tree
        std::vector<size_type> m_nodes_connects_size;
        std::vector<size_type> m_nodes_connects_ptr;
        std::vector<size_type> m_nodes_adjacency;
        std::vector<std::tuple<size_type /* node */,
                               size_type /* parent */,
                               data_type /* pass elevation */,
                               data_type /* parent pass elevation */>>
            m_reorder_stack;

        // TODO: make it an option (true for sink route fill sloped)
        bool m_keep_order = true;

        std::vector<size_type> m_pass_stack;
        std::vector<size_type> m_parent_basins;

        friend class ::BasinGraph_Test;
    };


    template <class FG>
    void basin_graph<FG>::compute_basins()
    {
        auto& basins = m_flow_graph_impl->m_basins;
        const auto& receivers = m_flow_graph_impl.receivers();
        const auto& dfs_indices = m_flow_graph_impl.dfs_indices();

        size_type current_basin;
        current_basin = -1;

        m_outlets.clear();

        for (const auto& idfs : dfs_indices)
        {
            if (idfs == receivers(idfs))
            {
                m_outlets.push_back(idfs);
                current_basin++;
            }

            basins(idfs) = current_basin;
        }

        assert(m_outlets.size() == current_basin + 1);
    }

    template <class FG>
    void basin_graph<FG>::update_routes(const data_array_type& elevation)
    {
        connect_basin(elevation);

        if (m_mst_method == mst_method::kruskal)
        {
            compute_tree_kruskal();
        }
        else
        {
            compute_tree_boruvka();
        }

        orient_edges();
    }

    template <class FG>
    void basin_graph<FG>::connect_basins(const data_array_type& elevation)
    {
        using neighbors_type = typename flow_graph_impl_type::grid_type::neighbors_type;

        auto nbasins = basins_count();

        auto& basins = m_flow_graph_impl->m_basins;
        const auto& receivers = m_flow_graph_impl.receivers();
        const auto& dfs_indices = m_flow_graph_impl.dfs_indices();

        auto& grid = m_flow_graph_impl.grid();
        const auto& status_at_nodes = grid.status_at_nodes();

        neighbors_type neighbors;

        size_type ibasin;
        size_type current_basin = -1;

        // assume iteration starts at a base level node (outer basin)!
        bool is_inner_basin = false;

        // reset
        m_root = -1;
        m_edges.clear();
        m_edges.reserve(4 * nbasins);

        m_edge_positions.resize(nbasins);
        std::fill(m_edge_positions.begin(), m_edge_positions.end(), -1);
        m_edge_positions_tmp.reserve(nbasins);
        m_edge_positions_tmp.clear();

        for (const auto idfs : dfs_indices)
        {
            const auto irec = receivers(idfs, 0);

            // any new basin visited
            if (irec == idfs)
            {
                ibasin = basins(idfs);
                is_inner_basin = status_at_nodes.flat(idfs) != node_status::fixed_value_boundary;

                if (!is_inner_basin)
                {
                    if (m_root == -1)
                    {
                        // assign root to an existing outer basin
                        m_root = ibasin;
                    }
                    else
                    {
                        m_edges.push_back(edge::make_edge(m_root, ibasin));
                    }
                }
            }

            // any node in an inner basin
            if (is_inner_basin)
            {
                const data_type ielev = elevation.flat(idfs);

                for (auto n : grid.neighbors(idfs, neighbors))
                {
                    const size_type nbasin = basins(n.idx);

                    // skip if neighbor node is in the same basin or in an
                    // already connected adjacent basin unless the latter is an
                    // outer basin
                    bool skip = ibasin >= nbasin;
                    node_status nstatus = status_at_nodes.flat(m_outlets[nbasin]);
                    bool is_inner_nbasin = nstatus != node_status::fixed_value_boundary;
                    if (skip && is_inner_nbasin)
                    {
                        continue;
                    }

                    const data_type pass_elevation = std::max(ielev, elevation.flat(n.idx));

                    // just jumped from one basin to another
                    // -> update current basin and reset its visited neighbor basins
                    if (current_basin != ibasin)
                    {
                        for (const auto& ivisited : m_edge_positions_tmp)
                        {
                            m_edge_positions[ivisited] = -1;
                        }
                        m_edge_positions_tmp.clear();
                        current_basin = ibasin;
                    }

                    // try getting the position (index) of the edge connecting
                    // with the adjacent basin:
                    // - if it is undefined (-1), add a new edge
                    // - if it is defined, update the edge if a pass of lower elevation is found
                    const size_type edge_idx = m_edge_positions[nbasin];

                    if (edge_idx == -1)
                    {
                        m_edge_positions[nbasin] = m_edges.size();
                        m_edge_positions_tmp.push_back(nbasin);

                        m_edges.push_back({ { ibasin, nbasin }, { idfs, n.idx }, pass_elevation });
                    }
                    else if (pass_elevation < m_edges[edge_idx].pass_elevation)
                    {
                        m_edges[edge_idx]
                            = edge{ { ibasin, nbasin }, { idfs, n.idx }, pass_elevation };
                    }
                }
            }
        }
    }

    template <class FG>
    void basin_graph<FG>::compute_tree_kruskal()
    {
        m_tree.reserve(basins_count() - 1);
        m_tree.clear();

        // sort edges indices by edge weight (pass elevation)
        m_edges_indices.resize(m_edges.size());
        std::iota(m_edges_indices.begin(), m_edges_indices.end(), 0);
        std::sort(m_edges_indices.begin(),
                  m_edges_indices.end(),
                  [&m_edges = m_edges](const size_type& i0, const size_type& i1)
                  { return m_edges[i0].pass_elevation < m_edges[i1].pass_elevation; });

        m_basins_uf.resize(basins_count());
        m_basins_uf.clear();

        for (size_type edge_idx : m_edges_indices)
        {
            size_type* link = m_edges[edge_idx].link;

            if (m_basins_uf.find(link[0]) != m_basins_uf.find(link[1]))
            {
                m_tree.push_back(edge_idx);
                m_basins_uf.merge(link[0], link[1]);
            }
        }
    }

    template <class FG>
    void basin_graph<FG>::compute_tree_boruvka()
    {
        const auto nbasins = basins_count();

        m_adjacency.clear();
        m_adjacency.resize(nbasins, { 0, 0 });
        m_low_degrees.reserve(nbasins);
        m_large_degrees.reserve(nbasins);

        m_edge_bucket.clear();
        m_edge_bucket.resize(nbasins, -1);

        // copy link basins
        m_link_basins.resize(m_edges.size());
        for (size_t i = 0; i < m_link_basins.size(); ++i)
        {
            m_link_basins[i][0] = m_edges[i].link[0];
            m_link_basins[i][1] = m_edges[i].link[1];
        }

        // first pass: create edge vector and compute adjacency size
        for (size_t lid = 0; lid < m_edges.size(); ++lid)
        {
            ++m_adjacency[m_link_basins[lid][0]].size;
            ++m_adjacency[m_link_basins[lid][1]].size;
        }

        // compute adjacency pointers
        m_adjacency[0].begin = 0;
        for (size_t nid = 1; nid < nbasins; ++nid)
        {
            m_adjacency[nid].begin = m_adjacency[nid - 1].begin + m_adjacency[nid - 1].size;
            m_adjacency[nid - 1].size = 0;
        }

        m_adjacency_list.resize(m_adjacency.back().begin + m_adjacency.back().size);
        m_adjacency.back().size = 0;

        for (size_t adj_data_i = 0; adj_data_i < m_adjacency_list.size(); ++adj_data_i)
            m_adjacency_list[adj_data_i].next = adj_data_i + 1;

        // second pass on edges: fill adjacency list
        for (size_t lid = 0; lid < m_edges.size(); ++lid)
        {
            auto& basins = m_link_basins[lid];

            m_adjacency_list[m_adjacency[basins[0]].begin + m_adjacency[basins[0]].size].link_id
                = lid;
            m_adjacency_list[m_adjacency[basins[1]].begin + m_adjacency[basins[1]].size].link_id
                = lid;

            ++m_adjacency[basins[0]].size;
            ++m_adjacency[basins[1]].size;
        }

        for (size_t nid = 0; nid < nbasins; ++nid)
        {
            check_capacity(m_low_degrees);
            check_capacity(m_large_degrees);
            // if degree is low enough
            if (m_adjacency[nid].size <= m_max_low_degree)
                m_low_degrees.push_back(nid);
            else
                m_large_degrees.push_back(nid);
        }

        m_perf_boruvka = 0;

        // compute the min span tree
        m_tree.reserve(nbasins - 1);
        m_tree.clear();

        while (m_low_degrees.size())
        {
            for (size_type nid : m_low_degrees)
            {
                // the node may have large degree after collapse
                if (m_adjacency[nid].size > m_max_low_degree)
                {
                    check_capacity(m_large_degrees);
                    m_large_degrees.push_back(nid);
                    continue;
                }

                // get the minimal weight edge that leaves that node
                size_type found_edge = -1;
                size_type node_B_id = -1;
                data_type found_edge_weight = std::numeric_limits<data_type>::max();

                size_type adjacency_data_ptr = m_adjacency[nid].begin;
                for (size_t step = 0; step < m_adjacency[nid].size; ++step)
                {
                    increase_perf_boruvka();

                    // find next adjacent edge in the list
                    size_type parsed_edge_id = m_adjacency_list[adjacency_data_ptr].link_id;
                    adjacency_data_ptr = m_adjacency_list[adjacency_data_ptr].next;

                    // check if the edge is valid (connected to a existing node)
                    // and if the weight is better than the previously found one
                    size_type opp_node = m_link_basins[parsed_edge_id][0];
                    if (opp_node == nid)
                        opp_node = m_link_basins[parsed_edge_id][1];

                    if (opp_node != nid && m_adjacency[opp_node].size > 0
                        && m_edges[parsed_edge_id].weight < found_edge_weight)
                    {
                        found_edge = parsed_edge_id;
                        found_edge_weight = m_edges[parsed_edge_id].weight;
                        node_B_id = opp_node;
                    }
                }

                if (found_edge == -1)
                    continue;  // TODO does it happens?

                // add edge to the tree
                check_capacity(m_tree);
                m_tree.push_back(found_edge);

                //  and collapse it toward opposite node

                // rename all A to B in adjacency of A
                adjacency_data_ptr = m_adjacency[nid].begin;
                for (size_t step = 0; step < m_adjacency[nid].size; ++step)
                {
                    increase_perf_boruvka();

                    // find next adjacent edge in the list
                    size_type edge_AC_id = m_adjacency_list[adjacency_data_ptr].link_id;

                    // TODO optimize that out?
                    if (step != m_adjacency[nid].size - 1)
                        adjacency_data_ptr = m_adjacency_list[adjacency_data_ptr].next;

                    // avoid self loop. A doesn't exist anymore, so edge AB
                    // will be discarded

                    if (m_link_basins[edge_AC_id][0] == nid)
                        m_link_basins[edge_AC_id][0] = node_B_id;
                    else
                        m_link_basins[edge_AC_id][1] = node_B_id;
                }

                // Append adjacency of B at the end of A
                m_adjacency_list[adjacency_data_ptr].next = m_adjacency[node_B_id].begin;

                // And collapse A into B
                m_adjacency[node_B_id].begin = m_adjacency[nid].begin;
                m_adjacency[node_B_id].size += m_adjacency[nid].size;

                // Remove the node from the graph
                m_adjacency[nid].size = 0;
            }

            m_low_degrees.clear();

            // Clean up graph (many edges are duplicates or self loops).
            int cur_large_degree = 0;
            for (size_type node_A_id : m_large_degrees)
            {
                // we will store all edges from A in the bucket, so that each edge
                // can appear only once
                m_edge_in_bucket.clear();
                size_type adjacency_data_ptr = m_adjacency[node_A_id].begin;

                for (size_t step = 0; step < m_adjacency[node_A_id].size; ++step)
                {
                    increase_perf_boruvka();

                    size_type edge_AB_id = m_adjacency_list[adjacency_data_ptr].link_id;
                    adjacency_data_ptr = m_adjacency_list[adjacency_data_ptr].next;

                    // find node B
                    size_type node_B_id = m_link_basins[edge_AB_id][0];
                    if (node_B_id == node_A_id)
                        node_B_id = m_link_basins[edge_AB_id][1];

                    if (m_adjacency[node_B_id].size > 0 && node_B_id != node_A_id)
                    {
                        // edge_bucket contain the edge_id connecting to opp_node_id
                        // or NodeId(-1)) if this is the first time we see it
                        size_type edge_AB_id_in_bucket = m_edge_bucket[node_B_id];

                        // first time we see
                        if (edge_AB_id_in_bucket == -1)
                        {
                            m_edge_bucket[node_B_id] = edge_AB_id;
                            m_edge_in_bucket.push_back(node_B_id);
                        }
                        else
                        {
                            // get weight of AB and of previously stored weight
                            data_type weight_in_bucket = m_edges[edge_AB_id_in_bucket].weight;
                            data_type weight_AB = m_edges[edge_AB_id].weight;

                            // if both weight are the same, we choose edge
                            // with min id
                            if (weight_in_bucket == weight_AB)
                                m_edge_bucket[node_B_id]
                                    = std::min(edge_AB_id_in_bucket, edge_AB_id);
                            else if (weight_AB < weight_in_bucket)
                                m_edge_bucket[node_B_id] = edge_AB_id;
                        }
                    }
                }

                // recompute connectivity of node A
                size_type cur_ptr = m_adjacency[node_A_id].begin;
                m_adjacency[node_A_id].size = m_edge_in_bucket.size();

                for (size_type node_B_id : m_edge_in_bucket)
                {
                    increase_perf_boruvka();

                    m_adjacency_list[cur_ptr].link_id = m_edge_bucket[node_B_id];
                    cur_ptr = m_adjacency_list[cur_ptr].next;

                    // clean occupency of edge_bucket for latter use
                    m_edge_bucket[node_B_id] = -1;
                }


                // update low degree information, if node A has low degree
                if (m_adjacency[node_A_id].size <= m_max_low_degree)
                {
                    check_capacity(m_low_degrees);
                    // add the node in low degree list
                    if (m_adjacency[node_A_id].size > 0)
                        m_low_degrees.push_back(node_A_id);
                }
                else
                    m_large_degrees[cur_large_degree++] = node_A_id;
            }
            m_large_degrees.resize(cur_large_degree);
            check_capacity(m_large_degrees);
        }
    }

    /*
     * Orient the edges of the basin graph (tree) in the upward (or counter)
     * flow direction.
     *
     * If needed, swap the link (i.e., basin node indices) and pass (i.e., grid
     * node indices) of the edges.
     */
    template <class FG>
    void basin_graph<FG>::orient_edges()
    {
        const auto nbasins = basins_count();

        // nodes connections
        m_nodes_connects_size.resize(nbasins);
        std::fill(m_nodes_connects_size.begin(), m_nodes_connects_size.end(), size_t(0));
        m_nodes_connects_ptr.resize(nbasins);

        // parse the edges to compute the number of edges per node
        for (size_type l_id : m_tree)
        {
            m_nodes_connects_size[m_edges[l_id].basins[0]] += 1;
            m_nodes_connects_size[m_edges[l_id].basins[1]] += 1;
        }

        // compute the id of first edge in adjacency table
        m_nodes_connects_ptr[0] = 0;
        for (size_t i = 1; i < nbasins; ++i)
        {
            m_nodes_connects_ptr[i] = (m_nodes_connects_ptr[i - 1] + m_nodes_connects_size[i - 1]);
            m_nodes_connects_size[i - 1] = 0;
        }

        // create the adjacency table
        m_nodes_adjacency.resize(m_nodes_connects_ptr.back() + m_nodes_connects_size.back());
        m_nodes_connects_size.back() = 0;

        // parse the edges to update the adjacency
        for (size_type l_id : m_tree)
        {
            size_type n0 = m_edges[l_id].basins[0];
            size_type n1 = m_edges[l_id].basins[1];
            m_nodes_adjacency[m_nodes_connects_ptr[n0] + m_nodes_connects_size[n0]] = l_id;
            m_nodes_adjacency[m_nodes_connects_ptr[n1] + m_nodes_connects_size[n1]] = l_id;
            m_nodes_connects_size[n0] += 1;
            m_nodes_connects_size[n1] += 1;
        }

        // depth-first parse of the tree, starting from basin0
        // stack of node, parent
        m_reorder_stack.reserve(nbasins);
        m_reorder_stack.clear();

        m_reorder_stack.push_back({ m_root,
                                    m_root,
                                    std::numeric_limits<data_type>::min(),
                                    std::numeric_limits<data_type>::min() });

        if (m_keep_order)
        {
            m_pass_stack.clear();
            m_parent_basins.resize(nbasins);
            m_parent_basins[m_root] = m_root;
        }

        while (m_reorder_stack.size())
        {
            size_type node, parent;
            data_type pass_elevation, parent_pass_elevation;
            std::tie(node, parent, pass_elevation, parent_pass_elevation) = m_reorder_stack.back();
            m_reorder_stack.pop_back();


            for (size_t i = m_nodes_connects_ptr[node];
                 i < m_nodes_connects_ptr[node] + m_nodes_connects_size[node];
                 ++i)
            {
                edge& edg = m_edges[m_nodes_adjacency[i]];

                // the edge comming from the parent node has already been updated.
                // in this case, the edge is (parent, node)
                if (edg.link[0] == parent && node != parent)
                {
                    if (m_keep_order)
                    {
                        // force children of base nodes to be parsed
                        if (pass_elevation <= parent_pass_elevation && edg.pass[0] != -1)
                            // the pass is bellow the water level of the parent basin
                            m_parent_basins[edg.link[1]] = m_parent_basins[edg.link[0]];
                        else
                        {
                            m_parent_basins[edg.link[1]] = edg.link[1];
                            m_pass_stack.push_back(m_nodes_adjacency[i]);
                        }
                    }
                }
                else
                {
                    // we want the edge to be (node, next), where next is upper in flow order
                    // we check if the first node of the edge is not "node"
                    if (node != edg.link[0])
                    {
                        std::swap(edg.link[0], edg.link[1]);
                        std::swap(edg.pass[0], edg.pass[1]);
                    }

                    m_reorder_stack.push_back({ edg.link[1],
                                                node,
                                                std::max(edg.pass_elevation, pass_elevation),
                                                pass_elevation });
                }
            }
        }
    }

    namespace detail
    {
        inline auto get_d8_distances_sep(double dx, double dy) -> std::array<double, 9>
        {
            std::array<double, 9> d8_dists;

            for (int k = 0; k < 9; ++k)
            {
                double d8_dx = dx * double(k % 3 - 1);
                double d8_dy = dy * double(k / 3 - 1);
                d8_dists[k] = std::sqrt(d8_dy * d8_dy + d8_dx * d8_dx);
            }

            return d8_dists;
        }

        template <class T>
        auto get_d8_distance_id(const T n1r, const T n1c, const T n2r, const T n2c)
        {
            T r = n1r - n2r + 1;
            T c = n1c - n2c + 1;
            return r + 3 * c;
        }

        template <class T>
        auto get_d8_distance_id(const T n1, const T n2, const T ncols)
        {
            return get_d8_distance_id(n1 / ncols, n1 % ncols, n2 / ncols, n2 % ncols);
        }

    }

    template <class FG>
    template <class Rcv_XT, class DistRcv_XT, class Elevation_XT>
    void basin_graph<FG>::update_pits_receivers(Rcv_XT& receivers,
                                                DistRcv_XT& dist2receivers,
                                                const Elevation_XT& elevation,
                                                double dx,
                                                double dy)
    {
        /* Update receivers of pit nodes (and possibly lowest pass nodes)
            based on basin connectivity.

            Distances to receivers are also updated. An infinite distance is
            arbitrarily assigned to pit nodes.

            A minimum spanning tree of the basin graph is used here. Edges of
            the graph are also assumed to be oriented in the inverse of flow direction.

        */

        const auto elev_shape = elevation.shape();
        const size_type nrows = (size_type) elev_shape[0];
        const size_type ncols = (size_type) elev_shape[1];

        const auto d8_distances = detail::get_d8_distances_sep(dx, dy);


        for (size_type l_id : m_tree)
        {
            edge& link = m_edges[l_id];

            // for readibility, hum...
#define OUTFLOW 0  // to
#define INFLOW 1   // from

            //    node_to = conn_nodes[i, 0]
            //  node_from = conn_nodes[i, 1]

            // skip open basins
            if (link.nodes[OUTFLOW] == -1)
                continue;

            size_type outlet_inflow = m_outlets[link.basins[INFLOW]];

            dist2receivers[outlet_inflow] = std::numeric_limits<double>::max();

            if (elevation[link.nodes[INFLOW]] < elevation[link.nodes[OUTFLOW]])
                receivers[outlet_inflow] = link.nodes[OUTFLOW];
            else
            {
                receivers[outlet_inflow] = link.nodes[INFLOW];
                receivers[link.nodes[INFLOW]] = link.nodes[OUTFLOW];

                dist2receivers(link.nodes[INFLOW]) = d8_distances[detail::get_d8_distance_id(
                    link.nodes[INFLOW], link.nodes[OUTFLOW], static_cast<size_type>(ncols))];
            }
        }
    }

    template <class FG>
    template <class Rcv_XT, class DistRcv_XT, class Elevation_XT>
    void basin_graph<FG>::update_pits_receivers_carve(Rcv_XT& receivers,
                                                      DistRcv_XT& dist2receivers,
                                                      const Elevation_XT& elevation,
                                                      double dx,
                                                      double dy)
    {
        const auto elev_shape = elevation.shape();
        const size_type nrows = (size_type) elev_shape[0];
        const size_type ncols = (size_type) elev_shape[1];

        const auto d8_distances = detail::get_d8_distances_sep(dx, dy);

        for (size_type l_id : m_tree)
        {
            edge& link = m_edges[l_id];
#define OUTFLOW 0  // to
#define INFLOW 1   // from

            // skip open basins
            if (link.nodes[OUTFLOW] == -1)
                continue;

            size_type outlet_inflow = m_outlets[link.basins[INFLOW]];
            size_type cur_node = link.nodes[INFLOW];
            size_type next_node = receivers(cur_node);
            data_type previous_dist = dist2receivers[cur_node];

            receivers(cur_node) = link.nodes[OUTFLOW];
            // std::cerr << "+ [" << cur_node << "]" << dist2receivers(cur_node);
            dist2receivers(cur_node) = d8_distances[detail::get_d8_distance_id(
                cur_node, link.nodes[OUTFLOW], static_cast<size_type>(ncols))];

            // std::cerr << "->" << dist2receivers(cur_node)<< std::endl;

            // std::cout << "Pass " << cur_node << " -> " << link.nodes[OUTFLOW] << std::endl;

            while (cur_node != outlet_inflow)
            {
                // std::cerr << "  [" << next_node << "]" << dist2receivers(next_node);
                std::swap(dist2receivers(next_node), previous_dist);
                // std::cerr << "->" << dist2receivers(cur_node)<< std::endl;


                size_type rcv_next_node = receivers(next_node);
                receivers(next_node) = cur_node;
                // std::cout << next_node << " -> " << cur_node << std::endl;
                cur_node = next_node;
                next_node = rcv_next_node;
            }
        }
    }

    template <class FG>
    template <class Rcv_XT, class DistRcv_XT, class Elevation_XT, class Basins_XT>
    void basin_graph<FG>::update_pits_receivers_sloped(Rcv_XT& receivers,
                                                       DistRcv_XT& dist2receivers,
                                                       const Elevation_XT& elevation,
                                                       const Basins_XT& basins,
                                                       double dx,
                                                       double dy)
    {
        const auto elev_shape = elevation.shape();
        const size_type nrows = (size_type) elev_shape[0];
        const size_type ncols = (size_type) elev_shape[1];

        const auto d8_distances = detail::get_d8_distances_sep(dx, dy);
        std::array<double, 9> d8_distances_inv;
        for (size_t i = 0; i < 9; ++i)
            d8_distances_inv[i] = 1.0 / d8_distances[i];

        std::queue<size_type> queue;

        enum class Tag : char
        {
            UnParsed = 0,
            InQueue = 1,
            WithRcv = 2
        };
        std::vector<Tag> tag(nrows * ncols, Tag::UnParsed);

        // parse in basin order
        for (const size_type pass : m_pass_stack)
        {
            const edge& l = m_edges[pass];
#define OUTFLOW 0  // to
#define INFLOW 1   // from

            if (l.nodes[OUTFLOW] == -1)
            {
                continue;
            }

            // receivers[l.nodes[INFLOW]] = l.nodes[OUTFLOW];
            // dist2receivers(l.nodes[INFLOW]) =
            // d8_distances[detail::get_d8_distance_id(l.nodes[INFLOW], l.nodes[OUTFLOW], ncols)];

            assert(tag[l.nodes[INFLOW]] == Tag::UnParsed);

            queue.push(l.nodes[INFLOW]);
            tag[l.nodes[OUTFLOW]] = Tag::WithRcv;
            tag[l.nodes[INFLOW]] = Tag::InQueue;
            const size_type parsed_basin = l.basins[INFLOW];
            assert(parsed_basin == m_parent_basins[parsed_basin]);

            auto outflow_coords = detail::coords(l.nodes[OUTFLOW], static_cast<size_type>(ncols));

            const data_type elev = l.pass_elevation;

            while (!queue.empty())
            {
                size_type node = queue.front();

                queue.pop();

                const auto coords = detail::coords(node, static_cast<size_type>(ncols));

                size_type rcv = -1;
                double rcv_cost = std::numeric_limits<double>::lowest();
                double cost_r = double(outflow_coords.first - coords.first);
                double cost_c = double(outflow_coords.second - coords.second);

                // parse neighbors
                for (int k = 1; k < 9; ++k)
                {
                    size_type rr = coords.first + fastscapelib::consts::d8_row_offsets[k];
                    size_type cc = coords.second + fastscapelib::consts::d8_col_offsets[k];

                    if (detail::in_bounds(elev_shape, rr, cc))
                    {
                        const size_type ineighbor = detail::index(rr, cc, ncols);

                        if ((ineighbor != l.nodes[OUTFLOW]
                             && m_parent_basins[basins(ineighbor)] != parsed_basin)
                            || elevation(ineighbor) > elev)
                            continue;

                        // neighbor is already parsed, in the same basin. Could be a receiver
                        if (tag[ineighbor] == Tag::WithRcv)
                        {
                            // cost is an angular distance to the outflow - node line.
                            double cost
                                = cost_r * double(fastscapelib::consts::d8_row_offsets[k])
                                  + cost_c * double(fastscapelib::consts::d8_col_offsets[k]);
                            cost *= d8_distances_inv[detail::get_d8_distance_id(
                                coords.first,
                                coords.second,
                                static_cast<size_type>(rr),
                                static_cast<size_type>(cc))];

                            if (cost > rcv_cost)
                            {
                                rcv = ineighbor;
                                rcv_cost = cost;
                            }
                        }

                        else if (tag[ineighbor] == Tag::UnParsed)
                        {
                            queue.push(ineighbor);
                            tag[ineighbor] = Tag::InQueue;
                        }
                    }
                }

                assert(rcv != -1);
                receivers(node) = rcv;
                dist2receivers(node) = d8_distances[detail::get_d8_distance_id(
                    node, rcv, static_cast<size_type>(ncols))];
                tag[node] = Tag::WithRcv;
            }
        }
    }

    template <class FG>
    template <mst_method algo,
              sink_route_method connect,
              class Basins_XT,
              class Rcv_XT,
              class DistRcv_XT,
              class Stack_XT,
              class Active_XT,
              class Elevation_XT>
    void basin_graph<FG>::update_receivers(Rcv_XT& receivers,
                                           DistRcv_XT& dist2receivers,
                                           const Basins_XT& basins,
                                           const Stack_XT& stack,
                                           const Active_XT& active_nodes,
                                           const Elevation_XT& elevation,
                                           data_type dx,
                                           data_type dy)
    {
        connect_basins(basins, receivers, stack, active_nodes, elevation);
        //        for(auto l : m_edges)
        //            std::cout << "[(" << l.basins[0] << ',' <<l.basins[1]<<")("
        //                      << l.nodes[0] << ',' <<l.nodes[1]<<") "<< l.weight << "] ";
        //        std::cout << std::endl;

        if (algo == mst_method::kruskal)
        {
            compute_tree_kruskal();
            //        for(auto t : _tree)
            //            std::cout << "(" << m_edges[t].basins[0] << ','
            //            <<m_edges[t].basins[1]<<")";
            //        std::cout << std::endl;
        }
        else
        {
            compute_tree_boruvka();
        }

        orient_edges();

        switch (connect)
        {
            case sink_route_method::basic:
                update_pits_receivers(receivers, dist2receivers, elevation, dx, dy);
                break;
            case sink_route_method::carve:
                update_pits_receivers_carve(receivers, dist2receivers, elevation, dx, dy);
                break;
            case sink_route_method::fill_sloped:
                update_pits_receivers_sloped(receivers, dist2receivers, elevation, basins, dx, dy);
                break;
            default:
                break;
        }
    }

}

#endif
