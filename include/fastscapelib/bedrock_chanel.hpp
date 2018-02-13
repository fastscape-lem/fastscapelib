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

#include "fastscapelib/basin_graph.hpp"

#include "fastscapelib/Profile.h"
#include <assert.h>
#include <xtensor/xio.hpp>

namespace fs = fastscapelib;


namespace fastscapelib
{

template<class Erosion_XT, class Elevation_XT, class Stack_XT,
         class Receivers_XT, class Dist2Receivers_XT,
         class Float_T = typename Erosion_XT::value_type>
void erode_spower(Erosion_XT& erosion, const Elevation_XT& elevation, const Stack_XT&stack,
                  const Receivers_XT& receivers, const Dist2Receivers_XT& dist2receivers,
                 Float_T area, Float_T k, Float_T m, Float_T n, Float_T dt, Float_T tolerance)
{

    const auto global_factor = k * dt * std::pow(area, m);

    for (const auto istack : stack)
    {
        const auto irec = receivers(istack);

        if (irec == istack)
        {
            // no erosion at basin outlets
            erosion(istack) = 0.0;
            continue;
        }

        const auto factor = global_factor / std::pow(dist2receivers(istack), n);

        const auto node_elevation = elevation(istack);
        const auto rcv_elevation = elevation(irec) - erosion(irec);

        // iterate: lower elevation until convergence
        auto elevation_k = node_elevation;
        auto elevation_prev = std::numeric_limits<Float_T>::max();

        while (std::abs(elevation_k - elevation_prev) > tolerance)
        {
            const auto slope = elevation_k - rcv_elevation;
            const auto diff = (elevation_k - node_elevation + factor * std::pow(slope, n)) /
                    (1. + factor * n * std::pow(slope, n - 1));
            elevation_k -= diff;
            elevation_prev = elevation_k;
        }

        erosion(istack) = node_elevation - elevation_k;
    }
}


}  // namespace fastscapelib
