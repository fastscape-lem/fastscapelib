#include "random_benchmark.hpp"

#include <xtensor/xrandom.hpp>

#include <chrono>

std::unordered_map<
    std::string,
    /*                     */ std::vector<std::pair<std::string, Benchmark::Func_Type>>>
    Benchmark::benchmarks;

static std::string
pptime(double time)
{
    std::stringstream ss;
    std::string unit = "ns";
    if (time >= 1000.0)
    {
        time /= 1000.0;
        unit = "us";

        if (time >= 1000.0)
        {
            time /= 1000.0;
            unit = "ms";

            if (time >= 1000.0)
            {
                time /= 1000.0;
                unit = "s";

                if (time >= 60.0)
                {
                    ss << std::floor(time / 60.0) << " min " << std::fmod(time, 60.0) << " sec";
                    return ss.str();
                }
            }
        }
    }

    ss << time << ' ' << unit;

    return ss.str();
}

void
Benchmark::do_run(std::string type,
                  size_t nrows,
                  size_t ncols,
                  int num_loops,
                  std::chrono::nanoseconds max_time,
                  bool pretty_print)
{
    for (auto F : benchmarks[type])
    {
        std::cout << F.first << std::endl;

        xt::random::seed(1);

        std::vector<int64_t> times;
        int64_t times_sum = 0;

        int loop_count;
        for (loop_count = 0; loop_count < num_loops; ++loop_count)
        {
            auto start = std::chrono::high_resolution_clock::now();

            F.second(nrows, ncols);

            auto stop = std::chrono::high_resolution_clock::now();

            times.push_back((stop - start).count());
            times_sum += times.back();

            if (times_sum > max_time.count() && loop_count >= 2)
                break;
        }

        if (loop_count < num_loops)
            ++loop_count;

        double times_avg = (double) times_sum / (double) times.size();
        double time_variance = 0;
        for (int64_t time : times)
            time_variance += ((double) time - times_avg) * ((double) time - times_avg);
        time_variance /= (double) (times.size() - 1);

        double time_sdev = std::sqrt(time_variance);

        if (pretty_print)
            std::cout << F.first << ": mean = " << pptime(times_avg)
                      << ", std-dev = " << pptime(time_sdev) << " (" << loop_count << ") loops"
                      << std::endl;
        else
            std::cout << F.first << " " << times_avg << " " << time_sdev << std::endl;
    }
}
