/**
 * @file
 * @brief helper class for graph oriented bassin filling
 * @author Guillaume Cordonnier
 *
 **/


#pragma once

#include <vector>
#include <limits>
#include <numeric>
#include <algorithm>
#include "fastscapelib/utils.hpp"
#include "fastscapelib/consts.hpp"
#include "fastscapelib/union_find.hpp"


namespace fs = fastscapelib;

namespace fastscapelib
{

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

    BasinGraph() {}

    template <class Basins_XT, class Stack_XT, class Rcv_XT>
    void compute_basins(Basins_XT& basins, const Stack_XT& stack,
                        const Rcv_XT& receivers);

    Basin_T basin_count() {return Basin_T(_outlets.size());}

    std::vector<Node_T>& outlets() {return _outlets;}


    template <class Basins_XT, class Rcv_XT, class Stack_XT,
              class Active_XT, class Elevation_XT>
    void connect_basins (const Basins_XT& basins, const Rcv_XT& receivers,
                         const Stack_XT& stack, const Active_XT& active_nodes,
                         const Elevation_XT& elevation);

    void compute_tree_kruskal();

    template<bool keep_order, class Elevation_XT>
    void reorder_tree(const Elevation_XT& elevation);

    template<class Rcv_XT, class DistRcv_XT, class Elevation_XT>
    void update_pits_receivers(Rcv_XT& receivers, DistRcv_XT& dist2receivers,
                               const Elevation_XT& elevation, double dx, double dy);
    void update_pits_receivers_continuous();
    void fill_sinks_flat();
    void fill_sinks_sloped();


    // Tests:
    bool is_links_eq(std::vector<Link_T>& oth)
    {
        if (oth.size() != _links.size())
            return false;
        for (size_t i = 0; i< _links.size(); ++i)
            if (!(_links[i] == oth[i]))
            {
                std::cout << '(' << _links[i].basins[0]<<'-'<<_links[i].basins[1] << "),(" << _links[i].nodes[0]<<'-'<<_links[i].nodes[1]  << ") w=" << _links[i].weight << std::endl;
                std::cout << '(' << oth[i].basins[0]<<'-'<<oth[i].basins[1] << "),(" << oth[i].nodes[0]<<'-'<<oth[i].nodes[1]  << ") w=" << oth[i].weight << std::endl;
                return false;
            }
        return true;
    }
    void print_links()
    {
        for (auto& l : _links)
            std::cout << '(' << l.basins[0]<<'-'<<l.basins[1] << "),(" << l.nodes[0]<<'-'<<l.nodes[1]  << ") w=" << l.weight << std::endl;
    }

protected:

    Basin_T add_basin(const Node_T& outlet)
    /*  */ {_outlets.push_back(outlet); return Basin_T(_outlets.size()-1);}



    //  add link to the link list
    void add_link(Link_T&& l) {_links.push_back(l);}

    // initialize the acceleration structures for update_link
    void init_update_link(Node_T nnodes)
    {
        // acceleration structures:
        conn_pos.resize(nnodes);
        std::fill(conn_pos.begin(), conn_pos.end(), (index_t)(-1));
        conn_pos_used.clear();

        conn_cur_basin = -1;
    }

    // change existing link only if it the weight is smaller
    // this is optimized, but should follow some rules:
    // basins nodes are sorted (basin[0] < basin[1])
    // and the function assumes that if basin[0] changes between two
    // invocations, then it will never be called with basin[0] again
    void update_link(Link_T&& l)
    {
        // check if first basin changed, in which case clean
        // the optimization structures
        const Basin_T& optim_basin = l.basins[0];
        const Basin_T& ineighbor_basin = l.basins[1];
        if (conn_cur_basin != optim_basin)
        {
            // clear used array
            for (auto& iused : conn_pos_used)
                iused = -1;
            conn_pos_used.clear();
            conn_cur_basin = optim_basin;
        }

        index_t& conn_idx = conn_pos[ineighbor_basin];
        if (conn_idx == -1)
        {
            conn_idx = _links.size(); /* reference */
            conn_pos_used.push_back(ineighbor_basin);
            add_link(std::move(l));
        }
        else if (l.weight < _links[conn_idx].weight)
            _links[conn_idx] = l;

    }



private:


    std::vector<Node_T> _outlets; // bottom nodes of basins
    std::vector<Link_T>  _links;
    std::vector<index_t> _tree; // indices of links

    Basin_T root;

    // acceleration for basin connections
    Basin_T conn_cur_basin;
    std::vector<index_t> conn_pos;
    std::vector<index_t> conn_pos_used;

    // kruskal
    std::vector<index_t> link_indices;
    UnionFind_T<Basin_T> basin_uf;

    //reorder tree
    std::vector<size_t> nodes_connects_size;
    std::vector<size_t> nodes_connects_ptr;
    std::vector<index_t> nodes_adjacency;
    std::vector<std::pair<Basin_T /*node*/, Basin_T /*parent*/>> reorder_stack;

    // reoder_tree, keep order
    std::vector<Link_T> passes;
    std::vector<Basin_T> basin_stack;

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

    const auto elev_shape = elevation.shape();
    const Node_T nrows = (Node_T) elev_shape[0];
    const Node_T ncols = (Node_T) elev_shape[1];
    const Node_T nnodes = nrows*ncols;

    auto flattened_elevation = detail::make_flattened(elevation);

    // root (sea basin)
    root = Basin_T(-1);
    Basin_T ibasin;

    init_update_link(nnodes);

    bool bactive = false;

    for (const auto& istack : stack)
    {
        const auto& irec = receivers(istack);

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
                    add_link(Link_T::outlink(root, ibasin));
            }
        }

        // any node on a inner basin
        if (bactive)
        {
            Node_T r, c; std::tie(r,c) = detail::coords(istack, ncols);

            for (int i = 1; i<5; ++i)
            {
                Node_T kr = r + (Node_T)consts::d4_row_offsets[i];
                Node_T kc = c + (Node_T)consts::d4_col_offsets[i];

                if (!detail::in_bounds(elev_shape, kr, kc))
                    continue;

                Node_T ineighbor = detail::index(kr, kc, ncols);
                const Basin_T& ineighbor_basin = basins(ineighbor);
                const Node_T&  ineighbor_outlet = _outlets[ineighbor];

                // skip same basin or already connected adjacent basin
                // don't skip adjacent basin if it's an open basin
                if (ibasin >= ineighbor_basin && active_nodes(ineighbor_outlet))
                    continue;

                update_link({{ibasin, ineighbor_basin},
                             {istack, ineighbor},
                             std::max(flattened_elevation(istack), flattened_elevation(ineighbor))});

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
              [&_links = _links](const Link_T& i0, const Link_T& i1) {return _links[i0].weight < _links[i1].weight;});

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
template<bool keep_order, class Elevation_XT>
void BasinGraph<Basin_T, Node_T, Elevation_T>::reorder_tree(const Elevation_XT& elevation)
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

    reorder_stack.push_back({root, root});

    // compile time if
    if (keep_order)
    {
        passes.resize(_outlets.size());
        passes[root] = Link_T::outlink(root, root);
        basin_stack.reserve(_outlets.size());
        basin_stack.clear();
        basin_stack.push_back(root);
    }

    while (reorder_stack.size())
    {
        Basin_T node, parent;
        std::tie(node, parent) = reorder_tree.back();
        reorder_tree.pop_back();


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
                    Link_T& pass = passes[node];
                    pass = link;

                    if(pass.nodes[0] != -1) //node is not a global min
                    {

                        Elevation_T pass_height =
                                std::max(
                                    elevation[pass.nodes[0]],
                                /**/elevation[pass.nodes[1]]);

                        // sill is below parent sill and parent is not
                        // a global min (otherwise flat starting areas are not handled))
                        if (pass_height <= passes[parent].weight
                                && passes[parent].nodes[0] != -1)
                        {
                            // in particular, basin[0] is set to parent
                            pass = passes[parent];
                        }
                        else
                        {
                            basin_stack.push_back(node);
                            // basin[0] is set to node, meaning that the
                            // found pass is at water level
                            pass.basins[0] = node;
                            pass.weight = pass_height;
                        }
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

                reorder_stack.push_back({link.basins[1], node});
            }
        }
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

            // update distance
            double delta_x = dx * double(link.nodes[INFLOW] % ncols - link.nodes[OUTFLOW] % ncols);
            double delta_y = dy * double(link.nodes[INFLOW] / ncols - link.nodes[OUTFLOW] / ncols);

            dist2receivers[link.nodes[INFLOW]] = std::sqrt(delta_x * delta_x + delta_y * delta_y);
        }
    }

}

}
