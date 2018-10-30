/**
 * @file
 * @brief helper class for graph oriented bassin filling
 * @author Guillaume Cordonnier
 *
 **/


#pragma once

#include <vector>
#include <array>
#include <queue>
#include <stack>
#include <limits>
#include <numeric>
#include <algorithm>
#include "fastscapelib/utils.hpp"
#include "fastscapelib/consts.hpp"
#include "fastscapelib/union_find.hpp"

#include "assert.h"
#include "fastscapelib/Profile.h"

class BasinGraph_Test;

namespace fastscapelib
{

enum class BasinAlgo {Kruskal, Boruvka};
enum class ConnectType {Simple, Carved, Sloped};

template <class Basin_T, class Node_T, class Weight_T>
struct Link
{
    Basin_T basins[2];
    Node_T  nodes[2];
    Weight_T weight;

    static Link outlink(const Basin_T& b0, const Basin_T& b1)
    {
        return Link{{b0, b1}, {Node_T(-1), Node_T(-1)},
            std::numeric_limits<Weight_T>::lowest()};
    }

    bool operator == (const Link& other)
    {
        return basins[0] == other.basins[0] && basins[1] == other.basins[1] &&
                nodes[0] == other.nodes[0]  && nodes[1] == other.nodes[1] &&
                weight == other.weight;
    }

};

template <class Basin_T, class Node_T, class Elevation_T = double>
class BasinGraph
{
public:

    using Link_T = Link<Basin_T, Node_T, Elevation_T>;

	BasinGraph() { perf_count_boruvka = -1; }



    template <class Basins_XT, class Stack_XT, class Rcv_XT>
    void compute_basins(Basins_XT& basins, const Stack_XT& stack,
                        const Rcv_XT& receivers);

    Basin_T basin_count() {return Basin_T(_outlets.size());}

    std::vector<Node_T>& outlets() {return _outlets;}



    template <BasinAlgo algo, ConnectType connect,
              class Basins_XT, class Rcv_XT, class DistRcv_XT,
              class Stack_XT, class Active_XT, class Elevation_XT>
    void update_receivers(Rcv_XT& receivers, DistRcv_XT& dist2receivers,
                          const Basins_XT& basins,
                          const Stack_XT& stack, const Active_XT& active_nodes,
                          const Elevation_XT& elevation, Elevation_T dx, Elevation_T dy);

	const std::vector<Link_T>&  getLinks() const { return _links; }
	const std::vector<index_t>& getTreeIndices() const { return _tree; }

	size_t getPerfBoruvka() const { return perf_count_boruvka; }

protected:


    Basin_T add_basin(const Node_T& outlet)
    /*  */ {_outlets.push_back(outlet); return Basin_T(_outlets.size()-1);}



    //  add link to the link list
    void add_link(Link_T&& l) {_links.push_back(l);}



    template <class Basins_XT, class Rcv_XT, class Stack_XT,
              class Active_XT, class Elevation_XT>
    void connect_basins (const Basins_XT& basins, const Rcv_XT& receivers,
                         const Stack_XT& stack, const Active_XT& active_nodes,
                         const Elevation_XT& elevation);

    void compute_tree_kruskal();
    template<int max_low_degree = 16> // 16 for d8, 8 for plannar graph
    void compute_tree_boruvka();

    template<bool keep_order>
    void reorder_tree();

    template<class Rcv_XT, class DistRcv_XT, class Elevation_XT>
    void update_pits_receivers(Rcv_XT& receivers, DistRcv_XT& dist2receivers,
                               const Elevation_XT& elevation, double dx, double dy);

    template<class Rcv_XT, class DistRcv_XT, class Elevation_XT>
    void update_pits_receivers_carve(Rcv_XT& receivers, DistRcv_XT& dist2receivers,
                                     const Elevation_XT& elevation, double dx, double dy);

    template<class Rcv_XT, class DistRcv_XT, class Elevation_XT, class Basins_XT>
    void update_pits_receivers_sloped(Rcv_XT& receivers, DistRcv_XT& dist2receivers,
                                      const Elevation_XT& elevation, const Basins_XT&, double dx, double dy);




private:

    std::vector<Node_T> _outlets; // bottom nodes of basins
    std::vector<Link_T>  _links;
    std::vector<index_t> _tree; // indices of links

    Basin_T root;

    // acceleration for basin connections
    std::vector<index_t> conn_pos;
    std::vector<index_t> conn_pos_used;

    // kruskal
    std::vector<index_t> link_indices;
    UnionFind_T<Basin_T> basin_uf;

    // boruvka
    std::vector<std::array<index_t, 2>> link_basins;
    struct Connect {index_t begin; size_t size;};
    std::vector<Connect> adjacency;
    struct EdgeParse {index_t link_id; index_t next;};
    std::vector<EdgeParse> adjacency_list;
    std::vector<index_t> low_degrees, large_degrees;
    std::vector<index_t> edge_bucket;
    std::vector<index_t> edge_in_bucket;

    //reorder tree
    std::vector<size_t> nodes_connects_size;
    std::vector<size_t> nodes_connects_ptr;
    std::vector<index_t> nodes_adjacency;
    std::vector<std::tuple<Basin_T /*node*/, Basin_T /*parent*/,
    /**/    Elevation_T /*pass height */, Elevation_T /* parent pass height */>> reorder_stack;

    // reoder_tree, keep order
    std::vector<index_t /*link id*/> pass_stack;
    std::vector<index_t> parent_basins;

    friend class ::BasinGraph_Test;

	size_t perf_count_boruvka;

};


template<class Basin_T, class Node_T, class Elevation_T>
template <class Basins_XT, class Stack_XT, class Rcv_XT>
void BasinGraph<Basin_T, Node_T, Elevation_T>::compute_basins(Basins_XT& basins,
                                                              const Stack_XT& stack,
                                                              const Rcv_XT& receivers)
{
    Basin_T cur_basin;
    cur_basin = -1;

    _outlets.clear();

    for(const auto& istack : stack)
    {
        if(istack == receivers(istack))
            cur_basin = add_basin(istack);

        basins(istack) = cur_basin;
    }
}

template<class Basin_T, class Node_T, class Elevation_T>
template <class Basins_XT, class Rcv_XT, class Stack_XT,
          class Active_XT, class Elevation_XT>
void BasinGraph<Basin_T, Node_T, Elevation_T>::connect_basins (const Basins_XT& basins, const Rcv_XT& receivers,
                                                               const Stack_XT& stack, const Active_XT& active_nodes,
                                                               const Elevation_XT& elevation)
{
    _links.clear();
    _links.reserve(4*_outlets.size());

    const auto elev_shape = elevation.shape();
    const Node_T ncols = (Node_T) elev_shape[1];

    // root (sea basin)
    root = Basin_T(-1);
    Basin_T ibasin;

    // acceleration structures:
    conn_pos.resize(_outlets.size());
    std::fill(conn_pos.begin(), conn_pos.end(), (index_t)(-1));
    conn_pos_used.reserve(_outlets.size());
    conn_pos_used.clear();

    Basin_T conn_cur_basin = -1;

    bool bactive = false;


    for (const auto istack : stack)
    {
        const auto irec = receivers(istack);

        // new basins
        if (irec == istack)
        {
            ibasin = basins(istack);
            bactive =  active_nodes(istack);

            if (!bactive)
            {
                if (root == Basin_T(-1))
                    root = ibasin;
                else
                {
                    add_link(Link_T::outlink(root, ibasin));
                }
            }
        }

        // any node on a inner basin
        if (bactive)
        {
            const auto rc = detail::coords(istack, ncols);
            const Node_T r = rc.first, c = rc.second;
            const Elevation_T elev = elevation(istack);

            for (int i = 1; i<9; ++i)
            {
                const Node_T kr = r + (Node_T)consts::d8_row_offsets[i];
                const Node_T kc = c + (Node_T)consts::d8_col_offsets[i];

                if (!detail::in_bounds(elev_shape, kr, kc))
                    continue;

                const Node_T ineighbor = detail::index(kr, kc, ncols);
                const Basin_T ineighbor_basin = basins(ineighbor);

                // skip same basin or already connected adjacent basin
                // don't skip adjacent basin if it is an open basin
                if (ibasin >= ineighbor_basin && active_nodes(_outlets[ineighbor_basin]))
                    continue;

                const Elevation_T weight = std::max(elev, elevation(ineighbor));


                if (conn_cur_basin != ibasin)
                {
                    // clear used array
                    for (const auto& iused : conn_pos_used)
                        conn_pos[iused] = -1;
                    conn_pos_used.clear();
                    conn_cur_basin = ibasin;
                }

                const index_t conn_idx = conn_pos[ineighbor_basin];
                if (conn_idx == -1)
                {
                    conn_pos[ineighbor_basin] = _links.size();
                    conn_pos_used.push_back(ineighbor_basin);

                    _links.push_back({{ibasin, ineighbor_basin},
                              {istack, ineighbor},
                              weight}
                            );
                }
                else if (weight < _links[conn_idx].weight)
                {
                    _links[conn_idx] = Link_T{{ibasin, ineighbor_basin},
                                        {istack, ineighbor},
                                        weight};
                }
            }

        }
    }
}

template<class Basin_T, class Node_T, class Elevation_T>
void BasinGraph<Basin_T, Node_T, Elevation_T>::compute_tree_kruskal()
{
    _tree.reserve(_outlets.size()-1);
    _tree.clear();

    // sort edges by indices
    link_indices.resize(_links.size());
    std::iota(link_indices.begin(), link_indices.end(), 0);
    std::sort(link_indices.begin(), link_indices.end(),
              [&_links = _links](const index_t& i0, const index_t& i1) {return _links[i0].weight < _links[i1].weight;});

    basin_uf.resize(_outlets.size());
    basin_uf.clear();

    for (index_t l_id : link_indices)
    {
        Basin_T* link_basins = _links[l_id].basins;

        if (basin_uf.find(link_basins[0]) != basin_uf.find(link_basins[1]))
        {
            _tree.push_back(l_id);
            basin_uf.merge(link_basins[0], link_basins[1]);
        }
    }
}

template<class Basin_T, class Node_T, class Elevation_T>
template<int max_low_degree> // 16 for d8, 8 for plannar graph
void BasinGraph<Basin_T, Node_T, Elevation_T>::compute_tree_boruvka()
{
	adjacency.clear();
	adjacency.resize(_outlets.size(), { 0,0 });
    low_degrees.reserve(_outlets.size());
    large_degrees.reserve(_outlets.size());

	edge_bucket.clear();
    edge_bucket.resize(_outlets.size(), -1);

    // copy link basins
    link_basins.resize(_links.size());
    for(size_t i = 0; i<link_basins.size(); ++i)
    {
        link_basins[i][0] = _links[i].basins[0];
        link_basins[i][1] = _links[i].basins[1];
    }

    // first pass: create edge vector and compute adjacency size
    for (size_t lid = 0; lid < _links.size(); ++lid)
    {
        ++adjacency[link_basins[lid][0]].size;
        ++adjacency[link_basins[lid][1]].size;
    }

    // compute adjacency pointers
    adjacency[0].begin = 0;
    for(size_t nid = 1; nid < _outlets.size(); ++nid)
    {
        adjacency[nid].begin = adjacency[nid-1].begin + adjacency[nid-1].size;
        adjacency[nid-1].size = 0;
    }

    adjacency_list.resize(adjacency.back().begin + adjacency.back().size);
    adjacency.back().size = 0;

    for (size_t adj_data_i = 0;  adj_data_i < adjacency_list.size(); ++adj_data_i)
        adjacency_list[adj_data_i].next = adj_data_i + 1;

    // second pass on edges: fill adjacency list
    for (size_t lid = 0; lid < _links.size(); ++lid)
    {
        auto& basins = link_basins[lid];

        adjacency_list[adjacency[basins[0]].begin + adjacency[basins[0]].size].link_id = lid;
        adjacency_list[adjacency[basins[1]].begin + adjacency[basins[1]].size].link_id = lid;

        ++adjacency[basins[0]].size;
        ++adjacency[basins[1]].size;
    }

#define CHECK_CAPACITY(a) assert(a.size() < a.capacity())

    for(size_t nid = 0; nid < _outlets.size(); ++nid)
    {
		CHECK_CAPACITY(low_degrees);
		CHECK_CAPACITY(large_degrees);
        // if degree is low enough
        if (adjacency[nid].size <= max_low_degree)
            low_degrees.push_back(nid);
        else
            large_degrees.push_back(nid);

    }

	perf_count_boruvka = 0;
#define PERF_INCREASE ++perf_count_boruvka


    // compute the min span tree
    _tree.reserve(_outlets.size()-1);
    _tree.clear();

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
            Elevation_T found_edge_weight = std::numeric_limits<Elevation_T>::max();

            index_t adjacency_data_ptr = adjacency[nid].begin;
            for(size_t step = 0; step < adjacency[nid].size; ++step)
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

                if (opp_node != nid && adjacency[opp_node].size > 0 &&
                        _links[parsed_edge_id].weight < found_edge_weight)
                {
                    found_edge = parsed_edge_id;
                    found_edge_weight = _links[parsed_edge_id].weight;
                    node_B_id = opp_node;
                }
            }

            if (found_edge == -1)
                continue; //TODO does it happens?

            // add edge to the tree
			CHECK_CAPACITY(_tree);
            _tree.push_back(found_edge);

            //  and collapse it toward opposite node

            // rename all A to B in adjacency of A
            adjacency_data_ptr = adjacency[nid].begin;
            for(size_t step = 0; step < adjacency[nid].size; ++step)
            {
				PERF_INCREASE;

                // find next adjacent edge in the list
                index_t edge_AC_id = adjacency_list[adjacency_data_ptr].link_id;

                // TODO optimize that out?
                if (step != adjacency[nid].size - 1)
                    adjacency_data_ptr = adjacency_list[adjacency_data_ptr].next;

                // avoid self loop. A doesn't exist anymore, so edge AB
                // will be discarded

                if(link_basins[edge_AC_id][0] == nid)
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
                    if(edge_AB_id_in_bucket == -1)
                    {
                        edge_bucket[node_B_id] = edge_AB_id;
                        edge_in_bucket.push_back(node_B_id);
                    }
                    else
                    {
                        // get weight of AB and of previously stored weight
                        Elevation_T weight_in_bucket = _links[edge_AB_id_in_bucket].weight;
                        Elevation_T weight_AB = _links[edge_AB_id].weight;

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

template<class Basin_T, class Node_T, class Elevation_T>
template<bool keep_order>
void BasinGraph<Basin_T, Node_T, Elevation_T>::reorder_tree()
{
    /*Orient the graph (tree) of basins so that the edges are directed in
        the inverse of the flow direction.

        If needed, swap values given for each edges (row) in `conn_basins`
        and `conn_nodes`.

    */

    // nodes connections
    nodes_connects_size.resize(_outlets.size());
    std::fill(nodes_connects_size.begin(), nodes_connects_size.end(), size_t(0));
    nodes_connects_ptr.resize(_outlets.size());

    // parse the edges to compute the number of edges per node
    for(index_t l_id: _tree)
    {
        nodes_connects_size[_links[l_id].basins[0]] += 1;
        nodes_connects_size[_links[l_id].basins[1]] += 1;
    }

    // compute the id of first edge in adjacency table
    nodes_connects_ptr[0] = 0;
    for (size_t i = 1; i<_outlets.size(); ++i)
    {
        nodes_connects_ptr[i] = (nodes_connects_ptr[i - 1] +
                nodes_connects_size[i - 1]);
        nodes_connects_size[i - 1] = 0;
    }

    // create the adjacency table
    nodes_adjacency.resize(nodes_connects_ptr.back() + nodes_connects_size.back());
    nodes_connects_size.back() = 0;

    // parse the edges to update the adjacency
    for (index_t l_id: _tree)
    {
        Basin_T n0 = _links[l_id].basins[0];
        Basin_T n1 = _links[l_id].basins[1];
        nodes_adjacency[nodes_connects_ptr[n0] + nodes_connects_size[n0]] = l_id;
        nodes_adjacency[nodes_connects_ptr[n1] + nodes_connects_size[n1]] = l_id;
        nodes_connects_size[n0] += 1;
        nodes_connects_size[n1] += 1;
    }

    // depth-first parse of the tree, starting from basin0
    // stack of node, parent
    reorder_stack.reserve(_outlets.size());
    reorder_stack.clear();

    reorder_stack.push_back({root, root,
                             std::numeric_limits<Elevation_T>::min(),
                             std::numeric_limits<Elevation_T>::min()});

    if (keep_order)
    {
        pass_stack.clear();
        parent_basins.resize(_outlets.size());
        parent_basins[root] = root;
    }

    while (reorder_stack.size())
    {
        Basin_T node, parent;
        Elevation_T pass_height, parent_pass_height;
        std::tie(node, parent, pass_height, parent_pass_height) = reorder_stack.back();
        reorder_stack.pop_back();


        for(size_t i = nodes_connects_ptr[node];
            i< nodes_connects_ptr[node] + nodes_connects_size[node];
            ++i)
        {
            Link_T& link = _links[nodes_adjacency[i]];

            // the edge comming from the parent node has already been updated.
            // in this case, the edge is (parent, node)
            if (link.basins[0] == parent && node != parent)
            {

                if (keep_order)
                {
					// force children of base nodes to be parsed
                    if(pass_height <= parent_pass_height && link.nodes[0] != -1)
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
                if(node != link.basins[0])
                {
                    std::swap(link.basins[0], link.basins[1]);
                    std::swap(link.nodes[0], link.nodes[1]);
                }

                reorder_stack.push_back({link.basins[1], node, std::max(link.weight, pass_height), pass_height});
            }
        }
    }
}

namespace detail
{
inline
auto get_d8_distances_sep(double dx, double dy) -> std::array<double, 9>
{
    std::array<double, 9> d8_dists;

    for(int k=0; k<9; ++k)
    {
        double d8_dx = dx * double(k % 3 -1);
        double d8_dy = dy * double(k / 3 -1);
        d8_dists[k] = std::sqrt(d8_dy*d8_dy + d8_dx*d8_dx);
    }

    return d8_dists;
}

template<class T>
auto get_d8_distance_id(const T n1r, const T n1c, const T n2r, const T n2c)
{
    T r = n1r - n2r +1;
    T c = n1c - n2c +1;
    return r + 3*c;
}

template<class T>
auto get_d8_distance_id(const T n1, const T n2, const T ncols)
{
    return get_d8_distance_id(n1 / ncols, n1 %ncols, n2 / ncols, n2%ncols);
}

}

template<class Basin_T, class Node_T, class Elevation_T>
template<class Rcv_XT, class DistRcv_XT, class Elevation_XT>
void BasinGraph<Basin_T, Node_T, Elevation_T>::update_pits_receivers(Rcv_XT& receivers, DistRcv_XT& dist2receivers,
                                                                     const Elevation_XT& elevation, double dx, double dy)
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


    for (index_t l_id : _tree)
    {
        Link_T& link = _links[l_id];

        // for readibility, hum...
#define OUTFLOW 0 // to
#define INFLOW  1 // from

        //    node_to = conn_nodes[i, 0]
        //  node_from = conn_nodes[i, 1]

        // skip open basins
        if (link.nodes[OUTFLOW] == -1)
            continue;

        Node_T outlet_inflow = _outlets[link.basins[INFLOW]];

        dist2receivers[outlet_inflow] = std::numeric_limits<double>::max();

        if(elevation[link.nodes[INFLOW]] < elevation[link.nodes[OUTFLOW]])
            receivers[outlet_inflow] = link.nodes[OUTFLOW];
        else
        {
            receivers[outlet_inflow] = link.nodes[INFLOW];
            receivers[link.nodes[INFLOW]] = link.nodes[OUTFLOW];

            dist2receivers(link.nodes[INFLOW]) = d8_distances[detail::get_d8_distance_id(link.nodes[INFLOW], link.nodes[OUTFLOW], ncols)];
        }
    }
}

template<class Basin_T, class Node_T, class Elevation_T>
template<class Rcv_XT, class DistRcv_XT, class Elevation_XT>
void BasinGraph<Basin_T, Node_T, Elevation_T>::update_pits_receivers_carve(Rcv_XT& receivers, DistRcv_XT& dist2receivers,
                                                                           const Elevation_XT& elevation, double dx, double dy)
{
    const auto elev_shape = elevation.shape();
    const index_t nrows = (index_t) elev_shape[0];
    const index_t ncols = (index_t) elev_shape[1];

    const auto d8_distances = detail::get_d8_distances_sep(dx, dy);

    for (index_t l_id : _tree)
    {
        Link_T& link = _links[l_id];
#define OUTFLOW 0 // to
#define INFLOW  1 // from

        // skip open basins
        if (link.nodes[OUTFLOW] == -1)
            continue;

        Node_T outlet_inflow = _outlets[link.basins[INFLOW]];
        Node_T cur_node = link.nodes[INFLOW];
        Node_T next_node = receivers(cur_node);
        Elevation_T previous_dist = dist2receivers[cur_node];

        receivers(cur_node) = link.nodes[OUTFLOW];
        //std::cerr << "+ [" << cur_node << "]" << dist2receivers(cur_node);
        dist2receivers(cur_node) = d8_distances[detail::get_d8_distance_id(cur_node, link.nodes[OUTFLOW], ncols)];
		
        //std::cerr << "->" << dist2receivers(cur_node)<< std::endl;

        //std::cout << "Pass " << cur_node << " -> " << link.nodes[OUTFLOW] << std::endl;

        while( cur_node != outlet_inflow)
        {
            //std::cerr << "  [" << next_node << "]" << dist2receivers(next_node);
            std::swap(dist2receivers(next_node), previous_dist);
            //std::cerr << "->" << dist2receivers(cur_node)<< std::endl;


            Node_T rcv_next_node = receivers(next_node);
            receivers(next_node) = cur_node;
            //std::cout << next_node << " -> " << cur_node << std::endl;
            cur_node = next_node;
            next_node = rcv_next_node;
        }
    }
}

template<class Basin_T, class Node_T, class Elevation_T>
template<class Rcv_XT, class DistRcv_XT, class Elevation_XT, class Basins_XT>
void BasinGraph<Basin_T, Node_T, Elevation_T>::update_pits_receivers_sloped(
    Rcv_XT& receivers,
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

    std::queue<Node_T> queue;

    enum class Tag : char { UnParsed = 0, InQueue = 1, WithRcv = 2};
    std::vector<Tag> tag (nrows * ncols, Tag::UnParsed);

    // parse in basin order
    for(const index_t pass : pass_stack)
    {
        const Link_T& l = _links[pass];
#define OUTFLOW 0 // to
#define INFLOW  1 // from

        if (l.nodes[OUTFLOW] == -1)
        {
            continue;
        }

        //receivers[l.nodes[INFLOW]] = l.nodes[OUTFLOW];
        //dist2receivers(l.nodes[INFLOW]) = d8_distances[detail::get_d8_distance_id(l.nodes[INFLOW], l.nodes[OUTFLOW], ncols)];

        assert(tag[l.nodes[INFLOW]] == Tag::UnParsed);

        queue.push(l.nodes[INFLOW]);
        tag[l.nodes[OUTFLOW]] = Tag::WithRcv;
        tag[l.nodes[INFLOW]] = Tag::InQueue;
        const Basin_T parsed_basin = l.basins[INFLOW];
        assert(parsed_basin == parent_basins[parsed_basin]);

        auto outflow_coords = detail::coords(l.nodes[OUTFLOW], ncols);

        const Elevation_T elev = l.weight;

        while (!queue.empty())
        {
            Node_T node = queue.front();

            queue.pop();

            const auto coords = detail::coords(node, ncols);

            Node_T rcv = -1;
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
                    const Node_T ineighbor = detail::index(rr, cc, ncols);

                    if ((ineighbor != l.nodes[OUTFLOW]
                         && parent_basins[basins(ineighbor)] != parsed_basin)
                        || elevation(ineighbor) > elev)
                        continue;

                    // neighbor is already parsed, in the same basin. Could be a receiver
                    if (tag[ineighbor] == Tag::WithRcv)
                    {
						// cost is an angular distance to the outflow - node line.
                        double cost = cost_r * double(fastscapelib::consts::d8_row_offsets[k]) + cost_c * double(fastscapelib::consts::d8_col_offsets[k]);
                        cost *= d8_distances_inv[detail::get_d8_distance_id(coords.first, coords.second, rr, cc)];

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
            dist2receivers(node) = d8_distances[detail::get_d8_distance_id(node, rcv, ncols)];
            tag[node] = Tag::WithRcv;
        }
    }
}

template<class Basin_T, class Node_T, class Elevation_T>
template <BasinAlgo algo, ConnectType connect,
          class Basins_XT, class Rcv_XT, class DistRcv_XT,
          class Stack_XT, class Active_XT, class Elevation_XT>
void BasinGraph<Basin_T, Node_T, Elevation_T>::update_receivers(
        Rcv_XT& receivers, DistRcv_XT& dist2receivers,
        const Basins_XT& basins,
        const Stack_XT& stack, const Active_XT& active_nodes,
        const Elevation_XT& elevation, Elevation_T dx, Elevation_T dy)
{
    {PROFILE(u0, "connect_basins");
        connect_basins(basins, receivers, stack, active_nodes, elevation);
//        for(auto l : _links)
//            std::cout << "[(" << l.basins[0] << ',' <<l.basins[1]<<")("
//                      << l.nodes[0] << ',' <<l.nodes[1]<<") "<< l.weight << "] ";
//        std::cout << std::endl;
    }
    if (algo == BasinAlgo::Kruskal)
    {PROFILE(u1, "compute_tree_kruskal");
        compute_tree_kruskal();
//        for(auto t : _tree)
//            std::cout << "(" << _links[t].basins[0] << ',' <<_links[t].basins[1]<<")";
//        std::cout << std::endl;
    }
    else
    {
        PROFILE(u1, "compute_tree_boruvka");
        compute_tree_boruvka();
    }
    {PROFILE(u2, "reorder_tree");
        reorder_tree<connect == ConnectType::Sloped>();
    }
    {PROFILE(u3, "update_pits_receivers");
        switch (connect) {
        case ConnectType::Simple:
            update_pits_receivers(receivers, dist2receivers,elevation, dx, dy);
            break;
        case ConnectType::Carved:
            update_pits_receivers_carve(receivers, dist2receivers,elevation, dx, dy);
            break;
        case ConnectType::Sloped:
            update_pits_receivers_sloped(receivers, dist2receivers,elevation,basins, dx, dy);
            break;
        default:
            break;
        }
        //std::cout << receivers << '\n' << dist2receivers << std::endl;

    }
}

}
