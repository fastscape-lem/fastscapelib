#pragma once

#include "benchmark.hpp"

#include "xtensor/xtensor.hpp"
#include <functional>

using FastscapeFunctionType = std::function<void(xt::xtensor<double, 2>&, xt::xtensor<bool, 2>&)>;

void fastscape_run(size_t, size_t, FastscapeFunctionType);


class RegisterFastscape
{
public:
    RegisterFastscape(std::string name, RandomFunctionType func)
        : _internal(name, "fastscape", std::bind(random_run, std::placeholders::_1, std::placeholders::_2, func)) {}

    Benchmark::Register _internal;
};
