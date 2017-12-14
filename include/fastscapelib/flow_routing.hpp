/**
 * @file
 * @brief Provides implementation of various efficient
 *        algorithms related to flow routing.
 */
#pragma once

#include <cmath>
#include <algorithm>
#include <limits>
#include <array>

#include "xtensor/xtensor.hpp"
#include "xtensor/xadapt.hpp"

#include "fastscapelib/utils.hpp"
#include "fastscapelib/consts.hpp"


namespace fs = fastscapelib;


namespace fastscapelib {

    namespace detail {

        inline auto get_d8_distances(double dx, double dy)
                                     -> std::array<double, 9> {
            std::array<double, 9> d8_dists;

            for(size_t k=0; k<9; ++k) {
                d8_dists[k] = std::sqrt(
                    std::pow(dy * fs::consts::d8_row_offsets[k], 2.0) +
                    std::pow(dx * fs::consts::d8_col_offsets[k], 2.0));
            }

            return d8_dists;
        }


        template<class A1, class A2, class A3>
        void add2stack(index_t& nstack,
                       A1& stack,
                       const A2& ndonors,
                       const A3& donors,
                       index_t inode) {
            for(unsigned short k=0; k<ndonors(inode); ++k) {
                index_t idonor = donors(inode, k);
                stack(nstack) = idonor;
                ++nstack;
                add2stack(nstack, stack, ndonors, donors, idonor);
            }
        }

    }

    template<class A1, class A2, class A3, class A4>
    void compute_receivers_d8(A1& receivers,
                              A2& dist2receivers,
                              const A3& elevation,
                              const A4& active_nodes,
                              double dx,
                              double dy) {
        using elev_t = typename A3::value_type;
        const auto d8_dists = detail::get_d8_distances(dx, dy);

        const auto elev_shape = elevation.shape();
        const index_t nrows = (index_t) elev_shape[0];
        const index_t ncols = (index_t) elev_shape[1];

        for(index_t r=0; r<nrows; ++r) {
            for(index_t c=0; c<ncols; ++c) {
                index_t inode = r * ncols + c;

                receivers(inode) = inode;
                dist2receivers(inode) = 0.;

                if(!active_nodes(r, c)) {
                    continue;
                }

                double slope_max = std::numeric_limits<double>::min();

                for(size_t k=1; k<=8; ++k) {
                    index_t kr = r + fs::consts::d8_row_offsets[k];
                    index_t kc = c + fs::consts::d8_col_offsets[k];

                    if(!fs::detail::in_bounds(elev_shape, kr, kc)) { continue; }

                    index_t ineighbor = kr * ncols + kc;
                    double slope = (elevation(r, c) - elevation(kr, kc)) / d8_dists[k];

                    if(slope > slope_max) {
                        slope_max = slope;
                        receivers(inode) = ineighbor;
                        dist2receivers(inode) = d8_dists[k];
                    }
                }
            }
        }

    }


    template<class A1, class A2, class A3>
    void compute_donors(A1& ndonors, A2& donors, const A3& receivers) {
        index_t nnodes = (index_t) receivers.size();

        std::fill(ndonors.begin(), ndonors.end(), 0);

        for(index_t inode=0; inode<nnodes; ++inode) {
            if(receivers(inode) != inode) {
                index_t irec = receivers(inode);
                donors(irec, ndonors(irec)) = inode;
                ++ndonors(irec);
            }
        }

    }


    template<class A1, class A2, class A3, class A4>
    void compute_stack(A1& stack,
                       const A2& ndonors,
                       const A3& donors,
                       const A4& receivers) {
        index_t nnodes = (index_t) receivers.size();
        index_t nstack = 0;

        for(index_t inode=0; inode<nnodes; ++inode) {
            if(receivers(inode) == inode) {
                stack(nstack) = inode;
                ++nstack;
                detail::add2stack(nstack, stack, ndonors, donors, inode);
            }
        }

    }


    template<class A1, class A2, class A3, class A4>
    index_t compute_basins(A1& basins,
                           A2& outlets,
                           const A3& stack,
                           const A4& receivers) {
        index_t ibasin = -1;

        for(auto&& istack : stack) {
            index_t irec = receivers(istack);

            if(irec == istack) {
                ++ibasin;
                outlets(ibasin) = istack;
            }

            basins(istack) = ibasin;
        }

        index_t nbasins = ibasin + 1;

        return nbasins;
    }


    template<class A1, class A2, class A3>
    index_t compute_pits(A1& pits,
                         const A2& outlets,
                         const A3& active_nodes,
                         index_t nbasins) {
        //TODO: works whether active_nodes is 1-d or 2-d but not safe!
        //      see xtensor issue #588
        index_t ipit = 0;

        for(index_t ibasin=0; ibasin<nbasins; ++ibasin) {
            index_t inode = outlets(ibasin);

            if(active_nodes(inode)) {
                pits(ipit) = inode;
                ++ipit;
            }
        }

        index_t npits = ipit;

        return npits;
    }


    template<class A1, class A2, class A3>
    void compute_drainage_area(A1& area, const A2& stack, const A3& receivers) {

        for(auto inode=stack.crbegin(); inode!=stack.crend(); ++inode) {
            if(receivers(*inode) != *inode) {
                area(receivers(*inode)) += area(*inode);
            }
        }
    }


    template<class A1, class A2, class A3>
    void compute_drainage_area(A1& area,
                               const A2& stack,
                               const A3& receivers,
                               double dx,
                               double dy) {
        std::fill(area.begin(), area.end(), dx * dy);

        //TODO: replace with safe/clean way to get flatten view in xtensor.
        //      (see xtensor issues #322 #324 #588).
        auto area_flat = xt::adapt(area.data(),
                                   std::array<size_t, 1>{ stack.size() });

        compute_drainage_area(area_flat, stack, receivers);
    }

}
