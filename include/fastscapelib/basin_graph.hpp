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

    static Link outlink(const Basin_T& b0, const Basin_T& b2)
    {
        return Link{{b0, b1}, {Node_T(-1), Node_T(-1)},
            std::numeric_limits<Weight_T>::lowest()};
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



protected:

    Basin_T add_basin(const Node_T& outlet)
    /*  */ {_basins.push_back(outlet); return Basin_T(_basins.size()-1);}
    Basin_T end_basin() {return Basin_T(_basins.size());}


    void add_link(Link_T&&);
    void update_link(Link_T&&);

private:


    std::vector<Basin_T> _outlets;
    std::vector<Link_T>  _links;

    Basin_T root;

    // acceleration for basin connections
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

    _outlets.clear();

    for(const auto& istack : stack)
    {
        if(istack == receivers(istack))
            cur_basin = basin_graph.add_basin(istack);

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

    detail::Flattened2D flattened_elevation(elevation);

    // root (sea basin)
    root = Basin_T(-1);
    Basin_T ibasin;

    // acceleration structures:
    conn_pos.resize(nnodes);
    conn_pos_used.clear();

    bool bactive = false;

    for (const auto& istack : stack)
    {
        const auto& irec = receivers(istack);

        // new basins
        if (irec == istack)
        {
            ibasin = basins(istack);
            bactive =  active_nodes(istack);

            // clear used array
            for (auto& iused : conn_pos_used)
                iused = -1;

            if (!bactive)
            {
                if (root == Basin_T(-1))
                    root = ibasin;
                else
                    add_link(Link::outlink(root, ibasin));
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

                Node_T ineighbor = detail::index(r, c, ncols);
                const Basin_T& ineighbor_basin = basins(ineighbor);
                const Node_T&  ineighbor_outlet = _outlets(ineighbor);

                // skip same basin or already connected adjacent basin
                // don't skip adjacent basin if it's an open basin
                if (ibasin >= ineighbor_basin && active_nodes(ineighbor_outlet))
                    continue;

                update_link({{ibasin, ineighbor_basin},
                             {istack, ineightbor},
                             std::max(flattened_elevation(istack), flattened_elevation(ineighbor))});

            }

        }
    }
}


}
