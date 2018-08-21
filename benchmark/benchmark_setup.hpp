#pragma once

#include <math.h>
#include <array>
#include <tuple>

#include "xtensor/xmath.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"
#include "xtensor/xview.hpp"


namespace benchmark_setup
{


/*
 * Run benchmarks for different grid sizes.
 *
 * Assumes a square grid, i.e., the total number of grid points is s^2.
 *
 * Use this function with benchmark macros, e.g.,
 *    BENCHMARK(...)->Apply(grid_sizes<benchmark::kMillisecond>)
 */
template<benchmark::TimeUnit time_unit>
static void grid_sizes(benchmark::internal::Benchmark* bench) {
    std::array<int, 5> sizes {256, 512, 1024, 2048, 4096};

    for (int s : sizes)
    {
        bench->Args({s})->Unit(time_unit);
    }
}


enum class SurfaceType {cone, cone_inv, cone_noise, flat_noise, gauss, custom};


/*
 * A class for generating one of the following synthetic topography on
 * a ``n`` x ``n`` square grid:
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
        nnodes_side = n;
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
                         + xt::random::rand<T>(shape) * 5. / nnodes_side);
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
    int nnodes_side;
    std::array<size_t, 2> shape;
    std::tuple<xt::xtensor<double, 2>, xt::xtensor<double, 2>> grid;
    xt::xtensor<T, 2> elevation;

    auto get_elevation_cone()
    {
        return std::sqrt(2.) - xt::sqrt(xt::pow(std::get<0>(grid), 2) +
                                        xt::pow(std::get<1>(grid), 2));
    }

};

/*
 * Set fixed boundary conditions for each of the 4 sides of the grid.
 */
template<class A>
void set_fixed_boundary_faces(A&& active_nodes)
{
    auto active_nodes_ = xt::view(active_nodes, xt::all(), xt::all());
    active_nodes_ = true;

    const auto shape = active_nodes.shape();
    const std::array<size_t, 2> rows_idx {0, shape[0] - 1};
    const std::array<size_t, 2> cols_idx {0, shape[1] - 1};

    auto row_bounds = xt::view(active_nodes, xt::keep(rows_idx), xt::all());
    row_bounds = false;

    auto col_bounds = xt::view(active_nodes, xt::all(), xt::keep(cols_idx));
    col_bounds = false;
}


/*
 * A base setup common to various benchmarks.
 */
template<SurfaceType surf_type, class T>
struct FastscapeSetupBase
{
    int nnodes;
    double dx = 100.;
    double dy = 100.;
    xt::xtensor<T, 2> elevation;
    xt::xtensor<bool, 2> active_nodes;
    xt::xtensor<index_t, 1> receivers;
    xt::xtensor<T, 1> dist2receivers;
    xt::xtensor<index_t, 1> ndonors;
    xt::xtensor<index_t, 2> donors;
    xt::xtensor<index_t, 1> stack;
    xt::xtensor<index_t, 1> basins;

    FastscapeSetupBase(int n)
    {
        auto topo = SyntheticTopography<surf_type, T>(n);

        elevation = topo.get_elevation();

        active_nodes = xt::ones<bool>(elevation.shape());
        set_fixed_boundary_faces(active_nodes);

        nnodes = n * n;

        receivers = xt::empty<index_t>({nnodes});
        dist2receivers = xt::empty<T>({nnodes});
        ndonors = xt::empty<index_t>({nnodes});
        donors = xt::empty<index_t>({nnodes, 8});
        stack = xt::empty<index_t>({nnodes});
        basins = xt::empty<index_t>({nnodes});
    }
};

}  // namespace benchmark_setup
