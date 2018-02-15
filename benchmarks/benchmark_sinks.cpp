#include "benchmark.hpp"
#include "random_benchmark.hpp"
#include "fastscape_benchmark.hpp"

#include "fastscapelib/sinks.hpp"
#include "fastscapelib/flow_routing.hpp"



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

void benchmark_fastscape_sinks(
	xt::xtensor<index_t, 1>&      stack,
	xt::xtensor<index_t, 1>&      receivers,
	xt::xtensor<double, 1>&       dist2receviers,
	const xt::xtensor<double, 2>& elevation,
	const xt::xtensor<bool, 2>&   active_nodes,
	double dx, double dy)
{

    //std::cout<< elevation << std::endl;

    xt::xtensor<double, 2> tmp_elevation = elevation;
    xt::xtensor<index_t, 1> ndonors(stack.shape());
    xt::xtensor<index_t, 2> donors(std::array<size_t, 2>{stack.shape()[0], 8});
	fs::fill_sinks_sloped(tmp_elevation);

    //std::cout<< tmp_elevation << std::endl;


	fs::compute_receivers_d8(receivers, dist2receviers, tmp_elevation, active_nodes, dx, dy);

    //std::cout<< receivers << std::endl;

	fs::compute_donors(ndonors, donors, receivers);

    //std::cout<< ndonors << std::endl;
    //std::cout<< donors << std::endl;


	fs::compute_stack(stack, ndonors, donors, receivers);



}

RegisterFastscape register_fastscape_sinks("Sinks", benchmark_fastscape_sinks);
