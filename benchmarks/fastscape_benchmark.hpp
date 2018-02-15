#pragma once

#include "benchmark.hpp"
#include "fastscapelib/utils.hpp"

#include "xtensor/xtensor.hpp"

#include <functional>

using FastscapeFunctionType = std::function<void(
	xt::xtensor<index_t, 1>&      /*stack*/, 
	xt::xtensor<index_t, 1>&      /*receivers*/,
	xt::xtensor<double, 1>&       /*dist2receviers*/,
	const xt::xtensor<double, 2>& /*elevation*/, 
	const xt::xtensor<bool, 2>&   /*active nodes  */,
	double dx, double dy)>;

void fastscape_run(size_t, size_t, FastscapeFunctionType);


class RegisterFastscape
{
public:
    RegisterFastscape(std::string name, FastscapeFunctionType func)
		: _internal(name, "fastscape", std::bind(fastscape_run, std::placeholders::_1, std::placeholders::_2, func)) {}

    Benchmark::Register _internal;
};
