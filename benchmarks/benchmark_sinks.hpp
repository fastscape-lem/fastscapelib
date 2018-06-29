#pragma once

void benchmark_fastscape_sinks(
	xt::xtensor<index_t, 1>&      stack,
	xt::xtensor<index_t, 1>&      receivers,
	xt::xtensor<double, 1>&       dist2receviers,
	const xt::xtensor<double, 2>& elevation,
	const xt::xtensor<bool, 2>&   active_nodes,
	double dx, double dy);