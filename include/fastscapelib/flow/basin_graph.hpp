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
                    mst_method mst_meth,
                    sink_route_method sink_route_meth)
            : m_flow_graph_impl(flow_graph_impl)
            , m_mst_method(mst_meth)
            , m_sink_route_method(sink_route_meth)
        {
            m_perf_boruvka = -1;

            // TODO: check shape of receivers (should be single flow)
        }

        size_type basins_count()
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
        const std::vector<index_t>& tree() const
        {
            return m_tree;
        }

        void compute_basins();

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

        template <int max_low_degree = 16>  // 16 for d8, 8 for plannar graph
        void compute_tree_boruvka();

        template <bool keep_order>
        void reorder_tree();

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
        std::vector<index_t> m_tree;  // indices of edges

        // root is used as a virtual basin "node" to which all outer basins are
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
        std::vector<index_t> edges_indices;
        UnionFind_T<size_type> basin_uf;

        // boruvka
        std::vector<std::array<index_t, 2>> link_basins;
        struct Connect
        {
            index_t begin;
            size_t size;
        };
        std::vector<Connect> adjacency;
        struct EdgeParse
        {
            index_t link_id;
            index_t next;
        };
        std::vector<EdgeParse> adjacency_list;
        std::vector<index_t> low_degrees, large_degrees;
        std::vector<index_t> edge_bucket;
        std::vector<index_t> edge_in_bucket;

        // reorder tree
        std::vector<size_t> nodes_connects_size;
        std::vector<size_t> nodes_connects_ptr;
        std::vector<index_t> nodes_adjacency;
        std::vector<std::tuple<size_type /*node*/,
                               size_type /*parent*/,
                               /**/ data_type /*pass height */,
                               data_type /* parent pass height */>>
            reorder_stack;

        // reoder_tree, keep order
        std::vector<index_t /*link id*/> pass_stack;
        std::vector<index_t> parent_basins;

        friend class ::BasinGraph_Test;

        size_t m_perf_boruvka;
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

        const auto elev_shape = elevation.shape();
        const size_type ncols = (size_type) elev_shape[1];

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
                        for (const auto& iused : m_edge_positions_tmp)
                        {
                            m_edge_positions[iused] = -1;
                        }
                        m_edge_positions_tmp.clear();
                        current_basin = ibasin;
                    }

                    // try getting the position (index) of the edge towards the adjacent basin:
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
        edges_indices.resize(m_edges.size());
        std::iota(edges_indices.begin(), edges_indices.end(), 0);
        std::sort(edges_indices.begin(),
                  edges_indices.end(),
                  [&m_edges = m_edges](const index_t& i0, const index_t& i1)
                  { return m_edges[i0].weight < m_edges[i1].weight; });

        basin_uf.resize(basins_count());
        basin_uf.clear();

        for (index_t edge_idx : edges_indices)
        {
            size_type* link = m_edges[edge_idx].link;

            if (basin_uf.find(link[0]) != basin_uf.find(link[1]))
            {
                m_tree.push_back(edge_idx);
                basin_uf.merge(link[0], link[1]);
            }
        }
    }

    template <class FG>
    template <int max_low_degree>  // 16 for d8, 8 for plannar graph
    void basin_graph<FG>::compute_tree_boruvka()
    {
        adjacency.clear();
        adjacency.resize(m_outlets.size(), { 0, 0 });
        low_degrees.reserve(m_outlets.size());
        large_degrees.reserve(m_outlets.size());

        edge_bucket.clear();
        edge_bucket.resize(m_outlets.size(), -1);

        // copy link basins
        link_basins.resize(m_edges.size());
        for (size_t i = 0; i < link_basins.size(); ++i)
        {
            link_basins[i][0] = m_edges[i].basins[0];
            link_basins[i][1] = m_edges[i].basins[1];
        }

        // first pass: create edge vector and compute adjacency size
        for (size_t lid = 0; lid < m_edges.size(); ++lid)
        {
            ++adjacency[link_basins[lid][0]].size;
            ++adjacency[link_basins[lid][1]].size;
        }

        // compute adjacency pointers
        adjacency[0].begin = 0;
        for (size_t nid = 1; nid < m_outlets.size(); ++nid)
        {
            adjacency[nid].begin = adjacency[nid - 1].begin + adjacency[nid - 1].size;
            adjacency[nid - 1].size = 0;
        }

        adjacency_list.resize(adjacency.back().begin + adjacency.back().size);
        adjacency.back().size = 0;

        for (size_t adj_data_i = 0; adj_data_i < adjacency_list.size(); ++adj_data_i)
            adjacency_list[adj_data_i].next = adj_data_i + 1;

        // second pass on edges: fill adjacency list
        for (size_t lid = 0; lid < m_edges.size(); ++lid)
        {
            auto& basins = link_basins[lid];

            adjacency_list[adjacency[basins[0]].begin + adjacency[basins[0]].size].link_id = lid;
            adjacency_list[adjacency[basins[1]].begin + adjacency[basins[1]].size].link_id = lid;

            ++adjacency[basins[0]].size;
            ++adjacency[basins[1]].size;
        }

#define CHECK_CAPACITY(a) assert(a.size() < a.capacity())

        for (size_t nid = 0; nid < m_outlets.size(); ++nid)
        {
            CHECK_CAPACITY(low_degrees);
            CHECK_CAPACITY(large_degrees);
            // if degree is low enough
            if (adjacency[nid].size <= max_low_degree)
                low_degrees.push_back(nid);
            else
                large_degrees.push_back(nid);
        }

        m_perf_boruvka = 0;
#define PERF_INCREASE ++m_perf_boruvka


        // compute the min span tree
        m_tree.reserve(m_outlets.size() - 1);
        m_tree.clear();

        while (low_degrees.size())
        {
            for (index_t nid : low_degrees)
            {
                // the node may have large degree after collapse
                if (adjacency[nid].size > max_low_degree)
                {
                    CHECK_CAPACITY(large_degrees);
                    large_degrees.push_back(nid);
                    continue;
                }

                // get the minimal weight edge that leaves that node
                index_t found_edge = -1;
                index_t node_B_id = -1;
                data_type found_edge_weight = std::numeric_limits<data_type>::max();

                index_t adjacency_data_ptr = adjacency[nid].begin;
                for (size_t step = 0; step < adjacency[nid].size; ++step)
                {
                    PERF_INCREASE;

                    // find next adjacent edge in the list
                    index_t parsed_edge_id = adjacency_list[adjacency_data_ptr].link_id;
                    adjacency_data_ptr = adjacency_list[adjacency_data_ptr].next;

                    // check if the edge is valid (connected to a existing node)
                    // and if the weight is better than the previously found one
                    index_t opp_node = link_basins[parsed_edge_id][0];
                    if (opp_node == nid)
                        opp_node = link_basins[parsed_edge_id][1];

                    if (opp_node != nid && adjacency[opp_node].size > 0
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
                CHECK_CAPACITY(m_tree);
                m_tree.push_back(found_edge);

                //  and collapse it toward opposite node

                // rename all A to B in adjacency of A
                adjacency_data_ptr = adjacency[nid].begin;
                for (size_t step = 0; step < adjacency[nid].size; ++step)
                {
                    PERF_INCREASE;

                    // find next adjacent edge in the list
                    index_t edge_AC_id = adjacency_list[adjacency_data_ptr].link_id;

                    // TODO optimize that out?
                    if (step != adjacency[nid].size - 1)
                        adjacency_data_ptr = adjacency_list[adjacency_data_ptr].next;

                    // avoid self loop. A doesn't exist anymore, so edge AB
                    // will be discarded

                    if (link_basins[edge_AC_id][0] == nid)
                        link_basins[edge_AC_id][0] = node_B_id;
                    else
                        link_basins[edge_AC_id][1] = node_B_id;
                }

                // Append adjacency of B at the end of A
                adjacency_list[adjacency_data_ptr].next = adjacency[node_B_id].begin;

                // And collapse A into B
                adjacency[node_B_id].begin = adjacency[nid].begin;
                adjacency[node_B_id].size += adjacency[nid].size;

                // Remove the node from the graph
                adjacency[nid].size = 0;
            }

            low_degrees.clear();

            // Clean up graph (many edges are duplicates or self loops).
            int cur_large_degree = 0;
            for (index_t node_A_id : large_degrees)
            {
                // we will store all edges from A in the bucket, so that each edge
                // can appear only once
                edge_in_bucket.clear();
                index_t adjacency_data_ptr = adjacency[node_A_id].begin;

                for (size_t step = 0; step < adjacency[node_A_id].size; ++step)
                {
                    PERF_INCREASE;

                    index_t edge_AB_id = adjacency_list[adjacency_data_ptr].link_id;
                    adjacency_data_ptr = adjacency_list[adjacency_data_ptr].next;

                    // find node B
                    index_t node_B_id = link_basins[edge_AB_id][0];
                    if (node_B_id == node_A_id)
                        node_B_id = link_basins[edge_AB_id][1];

                    if (adjacency[node_B_id].size > 0 && node_B_id != node_A_id)
                    {
                        // edge_bucket contain the edge_id connecting to opp_node_id
                        // or NodeId(-1)) if this is the first time we see it
                        index_t edge_AB_id_in_bucket = edge_bucket[node_B_id];

                        // first time we see
                        if (edge_AB_id_in_bucket == -1)
                        {
                            edge_bucket[node_B_id] = edge_AB_id;
                            edge_in_bucket.push_back(node_B_id);
                        }
                        else
                        {
                            // get weight of AB and of previously stored weight
                            data_type weight_in_bucket = m_edges[edge_AB_id_in_bucket].weight;
                            data_type weight_AB = m_edges[edge_AB_id].weight;

                            // if both weight are the same, we choose edge
                            // with min id
                            if (weight_in_bucket == weight_AB)
                                edge_bucket[node_B_id] = std::min(edge_AB_id_in_bucket, edge_AB_id);
                            else if (weight_AB < weight_in_bucket)
                                edge_bucket[node_B_id] = edge_AB_id;
                        }
                    }
                }

                // recompute connectivity of node A
                index_t cur_ptr = adjacency[node_A_id].begin;
                adjacency[node_A_id].size = edge_in_bucket.size();

                for (index_t node_B_id : edge_in_bucket)
                {
                    PERF_INCREASE;

                    adjacency_list[cur_ptr].link_id = edge_bucket[node_B_id];
                    cur_ptr = adjacency_list[cur_ptr].next;

                    // clean occupency of edge_bucket for latter use
                    edge_bucket[node_B_id] = -1;
                }


                // update low degree information, if node A has low degree
                if (adjacency[node_A_id].size <= max_low_degree)
                {
                    CHECK_CAPACITY(low_degrees);
                    // add the node in low degree list
                    if (adjacency[node_A_id].size > 0)
                        low_degrees.push_back(node_A_id);
                }
                else
                    large_degrees[cur_large_degree++] = node_A_id;
            }
            large_degrees.resize(cur_large_degree);
            CHECK_CAPACITY(large_degrees);
        }
    }

    template <class FG>
    template <bool keep_order>
    void basin_graph<FG>::reorder_tree()
    {
        /*Orient the graph (tree) of basins so that the edges are directed in
            the inverse of the flow direction.

            If needed, swap values given for each edges (row) in `conn_basins`
            and `conn_nodes`.

        */

        // nodes connections
        nodes_connects_size.resize(m_outlets.size());
        std::fill(nodes_connects_size.begin(), nodes_connects_size.end(), size_t(0));
        nodes_connects_ptr.resize(m_outlets.size());

        // parse the edges to compute the number of edges per node
        for (index_t l_id : m_tree)
        {
            nodes_connects_size[m_edges[l_id].basins[0]] += 1;
            nodes_connects_size[m_edges[l_id].basins[1]] += 1;
        }

        // compute the id of first edge in adjacency table
        nodes_connects_ptr[0] = 0;
        for (size_t i = 1; i < m_outlets.size(); ++i)
        {
            nodes_connects_ptr[i] = (nodes_connects_ptr[i - 1] + nodes_connects_size[i - 1]);
            nodes_connects_size[i - 1] = 0;
        }

        // create the adjacency table
        nodes_adjacency.resize(nodes_connects_ptr.back() + nodes_connects_size.back());
        nodes_connects_size.back() = 0;

        // parse the edges to update the adjacency
        for (index_t l_id : m_tree)
        {
            size_type n0 = m_edges[l_id].basins[0];
            size_type n1 = m_edges[l_id].basins[1];
            nodes_adjacency[nodes_connects_ptr[n0] + nodes_connects_size[n0]] = l_id;
            nodes_adjacency[nodes_connects_ptr[n1] + nodes_connects_size[n1]] = l_id;
            nodes_connects_size[n0] += 1;
            nodes_connects_size[n1] += 1;
        }

        // depth-first parse of the tree, starting from basin0
        // stack of node, parent
        reorder_stack.reserve(m_outlets.size());
        reorder_stack.clear();

        reorder_stack.push_back({ m_root,
                                  m_root,
                                  std::numeric_limits<data_type>::min(),
                                  std::numeric_limits<data_type>::min() });

        if (keep_order)
        {
            pass_stack.clear();
            parent_basins.resize(m_outlets.size());
            parent_basins[m_root] = m_root;
        }

        while (reorder_stack.size())
        {
            size_type node, parent;
            data_type pass_height, parent_pass_height;
            std::tie(node, parent, pass_height, parent_pass_height) = reorder_stack.back();
            reorder_stack.pop_back();


            for (size_t i = nodes_connects_ptr[node];
                 i < nodes_connects_ptr[node] + nodes_connects_size[node];
                 ++i)
            {
                edge& link = m_edges[nodes_adjacency[i]];

                // the edge comming from the parent node has already been updated.
                // in this case, the edge is (parent, node)
                if (link.basins[0] == parent && node != parent)
                {
                    if (keep_order)
                    {
                        // force children of base nodes to be parsed
                        if (pass_height <= parent_pass_height && link.nodes[0] != -1)
                            // the pass is bellow the water level of the parent basin
                            parent_basins[link.basins[1]] = parent_basins[link.basins[0]];
                        else
                        {
                            parent_basins[link.basins[1]] = link.basins[1];
                            pass_stack.push_back(nodes_adjacency[i]);
                        }
                    }
                }
                else
                {
                    // we want the edge to be (node, next), where next is upper in flow order
                    // we check if the first node of the edge is not "node"
                    if (node != link.basins[0])
                    {
                        std::swap(link.basins[0], link.basins[1]);
                        std::swap(link.nodes[0], link.nodes[1]);
                    }

                    reorder_stack.push_back({ link.basins[1],
                                              node,
                                              std::max(link.pass_elevation, pass_height),
                                              pass_height });
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
        const index_t nrows = (index_t) elev_shape[0];
        const index_t ncols = (index_t) elev_shape[1];

        const auto d8_distances = detail::get_d8_distances_sep(dx, dy);


        for (index_t l_id : m_tree)
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
        const index_t nrows = (index_t) elev_shape[0];
        const index_t ncols = (index_t) elev_shape[1];

        const auto d8_distances = detail::get_d8_distances_sep(dx, dy);

        for (index_t l_id : m_tree)
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
        const index_t nrows = (index_t) elev_shape[0];
        const index_t ncols = (index_t) elev_shape[1];

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
        for (const index_t pass : pass_stack)
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
            assert(parsed_basin == parent_basins[parsed_basin]);

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
                    index_t rr = coords.first + fastscapelib::consts::d8_row_offsets[k];
                    index_t cc = coords.second + fastscapelib::consts::d8_col_offsets[k];

                    if (detail::in_bounds(elev_shape, rr, cc))
                    {
                        const size_type ineighbor = detail::index(rr, cc, ncols);

                        if ((ineighbor != l.nodes[OUTFLOW]
                             && parent_basins[basins(ineighbor)] != parsed_basin)
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

        reorder_tree<connect == sink_route_method::fill_sloped>();

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
