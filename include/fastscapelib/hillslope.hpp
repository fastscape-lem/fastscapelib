/**
 * Functions to compute hillslope erosion.
 */
#pragma once

#include <array>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "xtensor/xbroadcast.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xnoalias.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xstrided_view.hpp"

#include "fastscapelib/utils.hpp"


namespace fastscapelib
{

namespace detail
{


/*
 * Solve tri-diagonal system of equations using Thomas' algorithm (TDMA).
 */
template<class L, class D, class U, class V>
auto solve_tridiagonal(L& lower,
                       D& diag,
                       U& upper,
                       V& vec)
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


/*
 * Get factors for linear diffusion ADI with a spatially uniform
 * diffusion coefficent (i.e., any scalar like float, double, etc.).
 */
template<class T>
auto get_adi_factors(T k_coef,
                     double dt,
                     double dx,
                     double dy,
                     size_t nrows,
                     size_t ncols,
                     typename std::enable_if_t<std::is_arithmetic<T>::value>* = 0)
{
    auto fr = xt::empty<double>({ncols});
    auto fc = xt::empty<double>({ncols});

    fr.fill(k_coef * 0.5 * dt / (dy * dy));
    fc.fill(k_coef * 0.5 * dt / (dx * dx));

    std::array<size_t, 3> f_shape {3, nrows, ncols};

    // TODO: don't use eval when xbroadcast will have strides
    //       (required for transpose views: see xtensor #917)
    auto&& factors_row = xt::eval(xt::broadcast(std::move(fr), f_shape));
    auto&& factors_col = xt::eval(xt::broadcast(std::move(fc), f_shape));

    return std::make_pair(factors_row, factors_col);
}


/*
 * Get factors for linear diffusion ADI with a spatially variable
 * diffusion coefficent (i.e., any xt::xexpression).
 */
template<class K>
auto get_adi_factors(K&& k_coef,
                     double dt,
                     double dx,
                     double dy,
                     size_t nrows,
                     size_t ncols,
                     typename std::enable_if_t<xt::is_xexpression<K>::value>* = 0)
{
    std::array<size_t, 3> f_shape {3, nrows, ncols};

    auto factors_row = xt::empty<double>(f_shape);
    auto factors_col = xt::empty<double>(f_shape);

    double fr = 0.25 * dt / (dy * dy);
    double fc = 0.25 * dt / (dx * dx);

    for (index_t r=1; r<nrows-1; ++r)
    {
        for (index_t c=1; c<ncols-1; ++c)
        {
            factors_row(0, r, c) = fr     * (k_coef(r-1, c) + k_coef(r, c));
            factors_row(1, r, c) = fr / 2 * (k_coef(r-1, c) + 2 * k_coef(r, c) + k_coef(r+1, c));
            factors_row(2, r, c) = fr     * (k_coef(r, c)   + k_coef(r+1, c));

            factors_col(0, r, c) = fc     * (k_coef(r, c-1) + k_coef(r, c));
            factors_col(1, r, c) = fc / 2 * (k_coef(r, c-1) + 2 * k_coef(r, c) + k_coef(r, c+1));
            factors_col(2, r, c) = fc     * (k_coef(r, c)   + k_coef(r, c+1));
        }
    }

    return std::make_pair(factors_row, factors_col);
}


/*
 * Solve linear diffusion for the row direction.
 */
template<class Ei, class Fr, class Fc>
auto solve_diffusion_adi_row(Ei&& elevation,
                             Fr&& factors_row,
                             Fc&& factors_col,
                             size_t nrows,
                             size_t ncols)
{
    xt::xtensor<double, 2> elevation_out = elevation;

    auto vec = xt::empty<double>({ncols});
    auto lower = xt::empty_like(vec);
    auto diag = xt::empty_like(vec);
    auto upper = xt::empty_like(vec);

    for (index_t r=1; r<nrows-1; ++r)
    {
        // TODO use xt::view with xbroadcast (check xtensor #1036 #917)
        xt::noalias(lower) = -1 * xt::strided_view(factors_col, {0, r, xt::all()});
        xt::noalias(diag) = 1 + 2 * xt::strided_view(factors_col, {1, r, xt::all()});
        xt::noalias(upper) = -1 * xt::strided_view(factors_col, {2, r, xt::all()});

        for (index_t c=1; c<ncols-1; ++c)
        {
            vec(c) = ((1 - 2 * factors_row(1, r, c)) * elevation(r, c) +
                      factors_row(0, r, c) * elevation(r-1, c) +
                      factors_row(2, r, c) * elevation(r+1, c));
        }

        // boundary conditions
        auto lower_bounds = xt::view(lower, xt::keep(0, -1));
        lower_bounds = 0;
        auto diag_bounds = xt::view(diag, xt::keep(0, -1));
        diag_bounds = 1;
        auto upper_bounds = xt::view(upper, xt::keep(0, -1));
        upper_bounds = 0;
        auto vec_bounds = xt::view(vec, xt::keep(0, -1));
        vec_bounds = xt::view(elevation, r, xt::keep(0, -1));

        auto elevation_out_r = xt::view(elevation_out, r, xt::all());

        elevation_out_r = solve_tridiagonal(lower, diag, upper, vec);
    }

    return elevation_out;
}


/*
 * Hillslope erosion by linear diffusion implementation.
 *
 * Note: this assumes row-major layout.
 */
template <class Er, class El, class K>
void erode_linear_diffusion_impl(Er&& erosion,
                                 El&& elevation,
                                 K&& k_coef,
                                 double dt,
                                 double dx,
                                 double dy)
{
    auto shape = elevation.shape();
    auto nrows = shape[0];
    auto ncols = shape[1];

    // TODO: optimize 0-d xexpression or xcalar using static_cast<double>?
    //       and/or assert k_coef.shape == elevation.shape?

    auto factors = get_adi_factors(std::forward<K>(k_coef), dt, dx, dy,
                                   nrows, ncols);

    // solve for rows
    auto elevation_tmp = solve_diffusion_adi_row(
        std::forward<El>(elevation),
        factors.first, factors.second,
        nrows, ncols);

    // solve for cols (i.e., transpose)
    auto elevation_next = solve_diffusion_adi_row(
        xt::transpose(elevation_tmp),
        xt::transpose(factors.second, {0, 2, 1}),
        xt::transpose(factors.first, {0, 2, 1}),
        ncols, nrows);

    erosion = elevation - xt::transpose(elevation_next);
}



}  // namespace detail

}  // namespace fastscapelib
