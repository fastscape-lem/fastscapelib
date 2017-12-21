#pragma once

#include <vector>
#include "xtensor/xadapt.hpp"
#include <initializer_list>


namespace fs = fastscapelib;

namespace fastscapelib
{

template <class Basin_T, class Node_T, class Weight_T>
struct Link
{
    Basin_T basins[2];
    Node_T  nodes[2];
    Weight_T weight;
};

template <class Basin_T, class Node_T, class Elevation_T = double>
class BasinGraph
{
public:

    using Link_T = Link<Basin_T, Node_T, Elevation_T>;

    BasinGraph() {}

    template <class Basins_XT, class Stack_XT, class Rcv_XT>
    void compute_basins(Basins_XT& basins,const Stack_XT& stack,const Rcv_XT& receivers);



protected:


    void clear() {_outlets.clear();}

    Basin_T add_basin(const Node_T& outlet) {_basins.push_back(outlet); return Basin_T(_basins.size()-1);}
    Basin_T end_basin() {return Basin_T(_basins.size());}


    void update_link(Link_T&&);


private:
    std::vector<Basin_T> _outlets;
    std::vector<Link_T>  _links;

};

template<class Basin_T, class Node_T, class Elevation_T>
template <class Basins_XT, class Stack_XT, class Rcv_XT>
void BasinGraph<Basin_T, Node_T, Elevation_T>::compute_basins(Basins_XT& basins,const Stack_XT& stack,const Rcv_XT& receivers)
{
    Basin_T cur_basin;

    basin_graph.clear();

    for(auto&& istack : stack)
    {
        if(istack == receivers(istack))
            cur_basin = basin_graph.add_basin(istack);

        basins(istack) = cur_basin;
    }
}


}
