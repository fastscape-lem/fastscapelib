#pragma once

#include <vector>
#include <limits>
#include "fastscapelib/utils.hpp"
#include "fastscapelib/consts.hpp"



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

    template <class Basins_XT, class Rcv_XT, class Stack_XT,
              class Active_XT, class Elevation_XT>
    void connect_basins (const Basins_XT& basins, const Rcv_XT& receivers,
                         const Stack_XT& stack, const Active_XT& active_nodes,
                         const Elevation_XT& elevation);

    Basin_T basin_count() {return Basin_T(_outlets.size());}

    std::vector<Basin_T>& outlets() {return _outlets;}


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


    std::vector<Basin_T> _outlets;
    std::vector<Link_T>  _links;

    Basin_T root;

    // acceleration for basin connections
    Basin_T conn_cur_basin;
    std::vector<index_t> conn_pos;
    std::vector<index_t> conn_pos_used;

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

}
