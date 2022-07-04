#pragma once

#include "xtensor/xtensor.hpp"
#include <chrono>

#include <unordered_map>

class Benchmark
{
public:
    Benchmark() = delete;

    template <class Duration>
    static void run(std::string type,
                    size_t nrows,
                    size_t ncols,
                    int num_loops,
                    Duration max_time,
                    bool pretty_print)
    {
        do_run(type,
               nrows,
               ncols,
               num_loops,
               std::chrono::duration_cast<std::chrono::nanoseconds>(max_time),
               pretty_print);
    }


    using Func_Type = std::function<void(size_t, size_t)>;

    class Register
    {
    public:
        Register(std::string name, std::string type, Func_Type f)
        {
            benchmarks[type].push_back({ name, f });
        }
    };

private:
    static void do_run(std::string type,
                       size_t,
                       size_t,
                       int num_loops,
                       std::chrono::nanoseconds max_time,
                       bool pretty_print);

private:
    static std::unordered_map<
        std::string,
        /*                     */ std::vector<std::pair<std::string, Func_Type>>>
        benchmarks;
};
