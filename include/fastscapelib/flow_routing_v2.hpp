#pragma once

#include "fastscapelib/utils.hpp"
#include "fastscapelib/consts.hpp"
#include "fastscapelib/union_find.hpp"

#include <vector>
#include <limits>
#include <assert.h>

namespace fastscapelib
{

namespace detail {
inline auto get_d8_distances_half(double dx, double dy) -> std::array<double, 4>
{
    std::array<double, 4> d8_dists;

    for(size_t k=0; k<4; ++k)
    {
        double d8_dx = dx * consts::d8_col_offsets_half[k];
        double d8_dy = dy * consts::d8_row_offsets_half[k];
        d8_dists[k] = std::sqrt(d8_dy*d8_dy + d8_dx*d8_dx);
    }

    return d8_dists;


}

template<class ElevT>
inline bool is_lower_edge(std::pair<ElevT, ElevT> e0, std::pair<ElevT, ElevT> e1)
{
    // lexicographic
    if (e0.first < e1.first)
        return true;
    else if (e0.first == e1.first)
        return e0.second < e1.second;
    else
        return false;
}

}


template<class A1, class A2, class A3, class A4>
void compute_receivers_passes(A1& receivers,
                          A2& passes,
                          const A3& elevation,
                          const A4& active_nodes,
                          double dx,
                          double dy)
{
    using Node_t = typename A1::value_type;
    using Elev_t = typename A3::value_type;
    const auto d8_dists = detail::get_d8_distances_half(dx, dy);

    const auto elev_shape = elevation.shape();
    const size_t nrows = (size_t) elev_shape[0];
    const size_t ncols = (size_t) elev_shape[1];

    // edges
    size_t max_edges = 4*nrows*ncols - 3*(nrows + ncols) + 2;

    std::vector<std::pair<Node_t, Node_t>> edge_nodes;
    std::vector<std::pair<Elev_t, Elev_t>> edge_weights;
    std::vector<index_t> graph_edges;
    std::vector<std::pair<Node_t, Node_t>> passes_list;

    edge_nodes.reserve(max_edges+1);
    edge_weights.reserve(max_edges+1);
    graph_edges.reserve(max_edges);
    passes_list.reserve(max_edges);


    //set edge 0 as a boudary condition: max weight
    edge_nodes.push_back({-1,-1});
    edge_weights.push_back({std::numeric_limits<Elev_t>::max(), 0.0});

    // build edge lists
    for(size_t r=0; r<nrows; ++r)
        for(size_t c=0; c<ncols; ++c)
        {
            index_t inode = r * ncols + c;

            Elev_t node_elev = elevation(inode);

            for(size_t k=0; k<4; ++k)
            {
                index_t kr = r + consts::d8_row_offsets_half[k];
                index_t kc = c + consts::d8_col_offsets_half[k];

                // TODO check is this can be optimized
                if(!detail::in_bounds(elev_shape, kr, kc))
                    continue;

                index_t ineighbor = kr * ncols + kc;

                Elev_t nb_elev = elevation(ineighbor);
                Elev_t slope = std::abs(nb_elev - node_elev)/d8_dists[k];

                graph_edges.push_back(edge_nodes.size());
                edge_nodes.push_back({inode, ineighbor});
                edge_weights.push_back({std::max(node_elev, nb_elev), -slope});
            }
        }

    // Union find for labbeling components.
    // Using it makes the algo O(n alpha(n)) instead of O(n),
    // but I suspect the alpha(n) to be always below the constant we would have with a graph traversal.
    // TODO check...
    UnionFind_T<Node_t> uf;
    uf.reserve(nrows*ncols);

    std::vector<index_t> graph_nodes;
    graph_nodes.reserve(nrows*ncols);

    // connect all inactive nodes to the first inactive
    // that will automatically discard any edges between inactive nodes
    Node_t first_inactive = -1;

    // find first inactive node
    for(first_inactive = 0;
        first_inactive < nrows*ncols && active_nodes(first_inactive);
        ++first_inactive)
    {
        graph_nodes.push_back(first_inactive);
        uf.push_back(first_inactive);
    }

    graph_nodes.push_back(first_inactive);

    // connect all inactive nodes together
    for (size_t i = first_inactive; i<nrows*ncols; ++i)
        if(!active_nodes(i))
            uf.push_back(first_inactive);
        else
        {
            uf.push_back(i);
            graph_nodes.push_back(i);
        }

    // internal data structures
    std::vector<index_t> lowest_nb_edge(nrows*ncols);
    std::vector<bool> used (nrows*ncols);

    std::vector<std::vector<index_t>> bucket_0(nrows*ncols);
    std::vector<index_t> in_bucket_0;
    in_bucket_0.reserve(nrows*ncols);
    std::vector<bool> bucket_1(nrows*ncols, false);

    // successive improvements of the graph
    while (graph_nodes.size() > 1)
    {
        // reset internal structures
        for(Node_t node : graph_nodes)
        {
            lowest_nb_edge[node] = 0;
            used[node] = false;
        }


        // find lowest edges adjacent to the nodes
        for(size_t eid : graph_edges)
        {
            Node_t u = uf.find(edge_nodes[eid].first);
            Node_t v = uf.find(edge_nodes[eid].second);

            // this only appens for ill formed initial graph,
            // i.e. where some vertices are already in the same connected componnent
            if(u == v)
                continue;

            index_t& lower_u = lowest_nb_edge[u];
            index_t& lower_v = lowest_nb_edge[v];

            const auto& edge_w = edge_weights[eid];

            if(detail::is_lower_edge(edge_w, edge_weights[lower_u]))
                lower_u = eid;
            if(detail::is_lower_edge(edge_w, edge_weights[lower_v]))
                lower_v = eid;
        }

        // merge connected components and update receiver list
        for(Node_t node : graph_nodes)
        {
            index_t eid = lowest_nb_edge[node];

            assert(eid !=0);

            // node is the graph representative of one of the node of the edge
            // we find back the two grid nodes that makes the edge
            index_t u = edge_nodes[eid].first;
            index_t v = edge_nodes[eid].second;

            uf.merge(u, v);

            // update receivers
            // TODO possible optimisation? the receivers should be all set after the first pass
            if(elevation(v) > elevation(u))
                passes_list.push_back({u, v});
            else
                receivers(u) = v;
        }

        // clean edges: sort each edge based on their lower representent
        graph_nodes.clear();
        in_bucket_0.clear();
        for(size_t eid : graph_edges)
        {
            Node_t repr = std::min(uf.find(edge_nodes[eid].first), uf.find(edge_nodes[eid].second));
            if(bucket_0[repr].empty())
            {
                in_bucket_0.push_back(repr);
                graph_nodes.push_back(repr);
                used[repr] = true;
            }
            bucket_0[repr].push_back(eid);
        }

        graph_edges.clear();

        for(index_t inb0 : in_bucket_0)
        {
            for(index_t eid : bucket_0[inb0])
            {
                Node_t repr = std::max(uf.find(edge_nodes[eid].first), uf.find(edge_nodes[eid].second));

                if (!bucket_1[repr])
                {
                    if(!used[repr])
                    {
                        used[repr] = true;
                        graph_nodes.push_back(repr);
                    }

                    bucket_1[repr] = true;
                    graph_edges.push_back(eid);
                }
            }
            for(index_t eid : bucket_0[inb0])
            {
                Node_t repr = std::max(uf.find(edge_nodes[eid].first), uf.find(edge_nodes[eid].second));
                bucket_1[repr] = false;
            }
            bucket_0[inb0].clear();
        }
    }

    // construct passes
    for(auto& p : passes)
        p = -1;
    for(auto& p : passes_list)
        // if the pass is ot a receiver
        // and the lowest node is active
        if (receivers[p.second] != p.first && active_nodes[p.first])
        {
            passes[p.first] = p.second;
            passes[p.second] = p.first;
        }

    //TODO check that there is no false passes near local min
}

void compute_stack_with_passes_list()
{

}

}
