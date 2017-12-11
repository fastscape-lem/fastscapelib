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
#include "xtensor/xview.hpp"

#include "fastscape/utils.hpp"
#include "fastscape/consts.hpp"


namespace fs = fastscape;


namespace fastscape {

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

        //TODO: shape assertions here? or do it (or already done) in callers...

        const auto elev_shape = elevation.shape();
        const index_t nrows = (index_t) elev_shape[0];
        const index_t ncols = (index_t) elev_shape[1];

        for(index_t r=0; r<nrows; ++r) {
            for(index_t c=0; c<ncols; ++c) {
                elev_t inode = r * ncols + c;

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

        for(index_t inode=0; inode<ndonors.size(); ++inode) {
            ndonors(inode) = (index_t) 0;
        }

        for(index_t inode=0; inode<receivers.size(); ++inode) {
            if(receivers(inode) != inode) {
                index_t irec = receivers(inode);
                donors(irec, ndonors(irec)) = inode;
                ++ndonors(irec);
            }
        }
    }

}
