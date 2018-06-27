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
		class Area_XT,
		class Float_T>
void erode_spower(Erosion_XT& erosion, const Elevation_XT& elevation, const Stack_XT&stack,
                  const Receivers_XT& receivers, const Dist2Receivers_XT& dist2receivers,
                 const Area_XT& area, Float_T k, Float_T m, Float_T n, Float_T dt, Float_T tolerance)
{

		using Scalar_T = std::common_type_t<Float_T, typename Erosion_XT::value_type, typename Elevation_XT::value_type, typename Area_XT::value_type>;

	std::vector<Scalar_T> water(stack.size());

    for (const auto istack : stack)
    {
        const auto irec = receivers(istack);

        if (irec == istack)
        {
            // no erosion at basin outlets
            erosion(istack) = 0.0;
			water[istack] = elevation(istack);
            continue;
        }

        const auto factor = k * dt * std::pow(area(istack), m) / std::pow(dist2receivers(istack), n);

        const auto node_elevation = elevation(istack);
        const auto rcv_elevation = elevation(irec) - erosion(irec);

		const auto s0 = rcv_elevation - node_elevation;

		if (s0 >= 0.0)
		{
			erosion(istack) = 0.0;
			water[istack] = std::max(node_elevation, water[irec]);
			continue;
		}

		// solve for s = elevation - rcv_elevation (proxy for slope)
		// equation becomes factor * s^n + s + s0 = 0
		// with s > 0

		// first guess. Always above the root
		auto s = -s0;

		for(;;)
		{
			const auto kxn = factor * std::pow(s, n);
			const auto knxnm1 = n * kxn / s;
			const auto F = kxn + s + s0;
			const auto dFds = 1 + knxnm1;

			// second order? To be tested and profiled
			// const auto d2Fds2 = (n - 1) * knxnm1 / slope;
			// const auto diff = 2 * F * dFds / (2 * dFds * dFds - F * d2Fds2)

			const auto diff = F / dFds; // diff > 0
			s -= diff;

			if (diff <= tolerance)
				break;
		} 

		auto elev = rcv_elevation + s;

		if (elev < water[irec])
			elev = std::min(node_elevation, water[irec]);

        erosion(istack) = node_elevation - elev;
		water[istack] = std::max(elev, water[irec]);
    }
}


}  // namespace fastscapelib
