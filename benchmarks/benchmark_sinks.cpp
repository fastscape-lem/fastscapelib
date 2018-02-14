#include "benchmark.hpp"
#include "random_benchmark.hpp"

#include "fastscapelib/sinks.hpp"



void benchmark_sinks_flat(xt::xtensor<double, 2>& elevation, xt::xtensor<bool, 2>&)
{
    fs::fill_sinks_flat(elevation);
}

RegisterRandom register_sinks_flat("Sinks Flat", benchmark_sinks_flat);


void benchmark_sinks_sloped(xt::xtensor<double, 2>& elevation, xt::xtensor<bool, 2>&)
{
    fs::fill_sinks_sloped(elevation);
}

RegisterRandom register_sinks_sloped("Sinks Sloped", benchmark_sinks_sloped);
