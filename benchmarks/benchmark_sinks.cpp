#include "benchmark.hpp"
#include "random_benchmark.hpp"
#include "fastscape_benchmark.hpp"
#include "benchmark_sinks.hpp"

#include "fastscapelib/sinks.hpp"
#include "fastscapelib/flow_routing.hpp"


namespace fs = fastscapelib;

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
	xt::xtensor<bm_index, 1>&      stack,
	xt::xtensor<bm_index, 1>&      receivers,
	xt::xtensor<bm_scalar, 1>&       dist2receviers,
	const xt::xtensor<bm_scalar, 2>& elevation,
	const xt::xtensor<bool, 2>&   active_nodes,
	bm_scalar dx, bm_scalar dy)
{

    //std::cout<< elevation << std::endl;

    xt::xtensor<bm_scalar, 2> tmp_elevation = elevation;
    xt::xtensor<bm_index, 1> ndonors(stack.shape());
    xt::xtensor<bm_index, 2> donors(std::array<size_t, 2>{stack.shape()[0], 8});
	fs::fill_sinks_sloped(tmp_elevation);

    //std::cout<< tmp_elevation << std::endl;


	fs::compute_receivers_d8(receivers, dist2receviers, tmp_elevation, active_nodes, dx, dy);

    //std::cout<< receivers << std::endl;

	fs::compute_donors(ndonors, donors, receivers);

    //std::cout<< ndonors << std::endl;
    //std::cout<< donors << std::endl;


	fs::compute_stack(stack, ndonors, donors, receivers);



}

//RegisterFastscape register_fastscape_sinks("Sinks", benchmark_fastscape_sinks);


#ifdef ENABLE_RICHDEM

#include "fastscapelib/richdem.hpp"

void benchmark_wei2018_flat(xt::xtensor<double, 2>& elevation, xt::xtensor<bool, 2>&)
{
	fs::fill_sinks_wei2018(elevation);
}

RegisterRandom register_sinks_wei2018("Wei2018 Flat", benchmark_wei2018_flat);


//WARNING: the following algorithm does not solve depression with a slight slope, as it would be excpected for fastscape.
// I still added the fastscape function for comparison with our methods on a single step.
// Although the results is not assumed to be correct, the timing should be consistent.

void benchmark_fastscape_wei2018(
	xt::xtensor<bm_index, 1>&      stack,
	xt::xtensor<bm_index, 1>&      receivers,
	xt::xtensor<bm_scalar, 1>&       dist2receviers,
	const xt::xtensor<bm_scalar, 2>& elevation,
	const xt::xtensor<bool, 2>&   active_nodes,
	bm_scalar dx, bm_scalar dy)
{

	xt::xtensor<bm_scalar, 2> tmp_elevation = elevation;
	xt::xtensor<bm_index, 1> ndonors(stack.shape());
	xt::xtensor<bm_index, 2> donors(std::array<size_t, 2>{stack.shape()[0], 8});

	
	fs::fill_sinks_wei2018(tmp_elevation);
	//fs::resolve_flats_sloped(tmp_elevation);


	fs::compute_receivers_d8(receivers, dist2receviers, tmp_elevation, active_nodes, dx, dy);


	fs::compute_donors(ndonors, donors, receivers);


	fs::compute_stack(stack, ndonors, donors, receivers);



}

void benchmark_fastscape_zhou2016(
	xt::xtensor<bm_index, 1>&      stack,
	xt::xtensor<bm_index, 1>&      receivers,
	xt::xtensor<bm_scalar, 1>&       dist2receviers,
	const xt::xtensor<bm_scalar, 2>& elevation,
	const xt::xtensor<bool, 2>&   active_nodes,
	bm_scalar dx, bm_scalar dy)
{

	xt::xtensor<bm_scalar, 2> tmp_elevation = elevation;
	xt::xtensor<bm_index, 1> ndonors(stack.shape());
	xt::xtensor<bm_index, 2> donors(std::array<size_t, 2>{stack.shape()[0], 8});


	fs::fill_sinks_zhou2016(tmp_elevation);
	//fs::resolve_flats_sloped(tmp_elevation);


	fs::compute_receivers_d8(receivers, dist2receviers, tmp_elevation, active_nodes, dx, dy);


	fs::compute_donors(ndonors, donors, receivers);


	fs::compute_stack(stack, ndonors, donors, receivers);



}

#endif