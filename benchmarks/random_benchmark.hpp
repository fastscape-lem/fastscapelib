#pragma once

#include "xtensor/xtensor.hpp"
#include <chrono>

class RandomBenchmark
{
public:

    RandomBenchmark() = delete;

    template<class Duration>
    static void run(size_t nrows, size_t ncols, int num_loops,
                    Duration max_time, bool pretty_print)
    {
        do_run(nrows, ncols, num_loops,
               std::chrono::duration_cast<std::chrono::nanoseconds>(max_time),
               pretty_print);
    }

    using Func_Type = void (*) (xt::xtensor<double, 2>&, xt::xtensor<bool, 2>&);

    class Register
    {
    public:
        Register(std::string name, Func_Type f) {benchmarks.push_back({name, f});}
    };

private:
    static void do_run(size_t, size_t, int num_loops,
                       std::chrono::nanoseconds max_time, bool pretty_print);

private:
    static std::vector<std::pair<std::string, Func_Type>> benchmarks;
};


