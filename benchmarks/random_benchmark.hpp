#ifndef RANDOM_BENCHMARK_H
#define RANDOM_BENCHMARK_H

#include "xtensor/xtensor.hpp"

class RandomBenchmark
{
public:

    RandomBenchmark() = delete;

    static void run(size_t, size_t, int num_loops, bool pretty_print);

    using Func_Type = void (*) (xt::xtensor<double, 2>&);

    class Register
    {
      public:
        Register(std::string name, Func_Type f) {benchmarks.push_back({name, f});}
    };


private:
    static std::vector<std::pair<std::string, Func_Type>> benchmarks;
    };

#endif // RANDOM_BENCHMARK_H
