#pragma once

#include "benchmark.hpp"

#include "xtensor/xtensor.hpp"
#include <functional>

using RandomFunctionType = std::function<void(xt::xtensor<double, 2>&, xt::xtensor<bool, 2>&)>;

void random_run(size_t, size_t, RandomFunctionType);


class RegisterRandom
{
public:
    RegisterRandom(std::string name, RandomFunctionType func)
        : _internal(name, "random", std::bind(random_run, std::placeholders::_1, std::placeholders::_2, func)) {}

    Benchmark::Register _internal;
};
