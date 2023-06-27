#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <tuple>

#include <benchmark/benchmark.h>

#include "xtensor/xmath.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xrandom.hpp"
#include "xtensor/xview.hpp"

#include "fastscapelib/grid/raster_grid.hpp"
#include "fastscapelib/utils/utils.hpp"


namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace bench_setup
    {

        /*
         * Run benchmarks for different grid sizes.
         *
         * Assumes a square grid, i.e., the total number of grid points is s^2.
         *
         * Use this function with benchmark macros, e.g.,
         *    BENCHMARK(...)->Apply(grid_sizes<benchmark::kMillisecond>)
         */
        template <benchmark::TimeUnit time_unit>
        static void grid_sizes(benchmark::internal::Benchmark* bench)
        {
            std::array<int, 5> sizes{ 256, 512, 1024, 2048, 4096 };

            for (int s : sizes)
            {
                bench->Args({ s })->Unit(time_unit);
            }
        }

        /*
         * Run benchmarks for different grid sizes.
         *
         * Only use 2 different sizes to keep it simple and readable.
         */
        template <benchmark::TimeUnit time_unit>
        static void small_grid_sizes(benchmark::internal::Benchmark* bench)
        {
            std::array<int, 2> sizes{ 512, 2048 };

            for (int s : sizes)
            {
                bench->Args({ s })->Unit(time_unit);
            }
        }

        enum class surface_type
        {
            cone,
            cone_inv,
            cone_noise,
            flat_noise,
            gauss,
            custom
        };

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
        template <surface_type S, class T>
        class synthetic_topography_2d
        {
        public:
            using grid_type = fs::raster_grid;
            using elevation_xt_type = xt::xtensor<T, 2>;

            int seed = 0;

            synthetic_topography_2d(int n);

            grid_type grid();
            elevation_xt_type elevation();

        private:
            using shape_type = typename grid_type::shape_type;

            int m_n;
            shape_type m_shape;
            std::tuple<xt::xtensor<double, 2>, xt::xtensor<double, 2>> m_mg;

            auto get_elevation_cone();
        };


        template <surface_type S, class T>
        synthetic_topography_2d<S, T>::synthetic_topography_2d(int n)
            : m_n(n)
            , m_shape({ static_cast<std::size_t>(n), static_cast<std::size_t>(n) })
        {
            m_mg = xt::meshgrid(xt::linspace<double>(-1, 1, n), xt::linspace<double>(-1, 1, n));
        }

        template <surface_type S, class T>
        auto synthetic_topography_2d<S, T>::grid() -> grid_type
        {
            fs::node_status fv = fs::node_status::fixed_value;

            grid_type grid = fs::raster_grid::from_length(m_shape, { 2.0, 2.0 }, fv);

            return grid;
        }

        template <surface_type S, class T>
        auto synthetic_topography_2d<S, T>::elevation() -> elevation_xt_type
        {
            elevation_xt_type elevation;

            xt::random::seed(seed);

            switch (S)
            {
                case surface_type::cone:
                    elevation = get_elevation_cone();
                    break;

                case surface_type::cone_inv:
                    elevation = -get_elevation_cone();
                    break;

                case surface_type::cone_noise:
                    elevation = (get_elevation_cone() + xt::random::rand<T>(m_shape) * 5. / m_n);
                    break;

                case surface_type::flat_noise:
                    elevation = xt::random::rand<T>(m_shape);
                    break;

                case surface_type::gauss:
                    elevation = xt::exp(
                        -(xt::pow(std::get<0>(m_mg), 2) / 2. + xt::pow(std::get<1>(m_mg), 2) / 2.));
                    break;

                default:
                    elevation = xt::zeros<T>(m_shape);
                    break;
            }

            return elevation;
        }

        template <surface_type S, class T>
        auto synthetic_topography_2d<S, T>::get_elevation_cone()
        {
            return std::sqrt(2.)
                   - xt::sqrt(xt::pow(std::get<0>(m_mg), 2) + xt::pow(std::get<1>(m_mg), 2));
        }

        /*
         * Set fixed boundary conditions for each of the 4 sides of the grid.
         */
        template <class A>
        void set_fixed_boundary_faces(A&& active_nodes)
        {
            auto active_nodes_ = xt::view(active_nodes, xt::all(), xt::all());
            active_nodes_ = true;

            auto row_bounds = xt::view(active_nodes, xt::keep(0, -1), xt::all());
            row_bounds = false;
            auto col_bounds = xt::view(active_nodes, xt::all(), xt::keep(0, -1));
            col_bounds = false;
        }

        /*
         * A base setup common to various benchmarks.
         */
        template <surface_type surf_type, class T>
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
                auto topo = synthetic_topography_2d<surf_type, T>(n);

                elevation = topo.elevation();

                active_nodes = xt::ones<bool>(elevation.shape());
                set_fixed_boundary_faces(active_nodes);

                nnodes = n * n;

                receivers = xt::empty<index_t>({ nnodes });
                dist2receivers = xt::empty<T>({ nnodes });
                ndonors = xt::empty<index_t>({ nnodes });
                donors = xt::empty<index_t>({ nnodes, 8 });
                stack = xt::empty<index_t>({ nnodes });
                basins = xt::empty<index_t>({ nnodes });
            }
        };

    }
}
