/**
 * Functions to compute hillslope erosion.
 */
#pragma once

#include <stdexcept>

#include "xtensor/xbuilder.hpp"

#include "fastscapelib/utils.hpp"


namespace fastscapelib
{

namespace detail
{


/*
 * Solve tri-diagonal system of equations using Thomas' algorithm (TDMA).
 */
template<class L, class D, class U, class V>
auto solve_tridiagonal(L&& lower,
                       D&& diag,
                       U&& upper,
                       V&& vec)
{
    auto n = vec.size();

    auto result = xt::empty_like(vec);
    auto gam = xt::empty_like(vec);

    if (diag(0) == 0)
    {
        throw std::runtime_error("division by zero while solving tri-diagonal system");
    }

    auto bet = diag(0);
    result(0) = vec(0) / bet;

    for (index_t i=1; i<n; ++i)
    {
        gam(i) = upper(i-1) / bet;
        bet = diag(i) - lower(i) * gam(i);

        if (bet == 0)
        {
            throw std::runtime_error("division by zero while solving tri-diagonal system");
        }

        result(i) = (vec(i) - lower(i) * result(i-1)) / bet;
    }

    for (index_t i=n-2; i>-1; --i)
    {
        result(i) -= gam(i+1) * result(i+1);
    }

    return result;
}


}  // namespace detail

}  // namespace fastscapelib
