#pragma once

#include <math.h>
#include <array>
#include <tuple>

#include "xtensor/xmath.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"


namespace fixtures
{


enum class SurfaceType {cone, cone_inv, cone_noise, flat_noise, gauss};


/*
 * A class for generating one of the following synthetic topography on
 * a n x n square grid:
 *
 * cone
 *   Cone surface (smooth, regular slope, no depression)
 * cone_inv
 *   Inverted cone surface (a single, big depression)
 * cone_noise
 *   Cone surface with random perturbations so that the surface
 *   has many depressions of different sizes
 * flat_noise
 *   Nearly flat surface with random perturbations
 *   (many small depressions)
 * gauss
 *   Gaussian surface (smooth, no depression)
 */
template<SurfaceType surf_type, class T>
class SyntheticTopography
{
public:
    int seed = 0;

    SyntheticTopography(int n)
    {
        npoints_side = n;
        auto n_ = static_cast<size_t>(n);
        shape = {n_, n_};

        grid = xt::meshgrid(xt::linspace<double>(-1, 1, n),
                            xt::linspace<double>(-1, 1, n));
    }

    xt::xtensor<T, 2> get_elevation()
    {
        xt::random::seed(seed);

        switch (surf_type) {
        case SurfaceType::cone:
            elevation = get_elevation_cone();
            break;

        case SurfaceType::cone_inv:
            elevation = -get_elevation_cone();
            break;

        case SurfaceType::cone_noise:
            elevation = (get_elevation_cone()
                         + xt::random::rand<T>(shape) * 5. / npoints_side);
            break;

        case SurfaceType::flat_noise:
            elevation = xt::random::rand<T>(shape);
            break;

        case SurfaceType::gauss:
            elevation = xt::exp(-(xt::pow(std::get<0>(grid), 2) / 2.
                                  + xt::pow(std::get<1>(grid), 2) / 2.));
            break;

        default:
            elevation = xt::zeros<T>(shape);
            break;
        }

        return elevation;
    }

private:
    int npoints_side;
    std::array<size_t, 2> shape;
    std::tuple<xt::xtensor<double, 2>, xt::xtensor<double, 2>> grid;
    xt::xtensor<T, 2> elevation;

    auto get_elevation_cone()
    {
        return std::sqrt(2.) - xt::sqrt(xt::pow(std::get<0>(grid), 2) +
                                        xt::pow(std::get<1>(grid), 2));
    }

};

}  // namespace fixtures
