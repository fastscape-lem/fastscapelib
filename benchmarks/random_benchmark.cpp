#include "random_benchmark.hpp"

#include <xtensor/xrandom.hpp>

#include <chrono>

std::vector<std::pair<std::string, RandomBenchmark::Func_Type>>
RandomBenchmark::benchmarks = std::vector<std::pair<std::string, RandomBenchmark::Func_Type>>();

static std::string pptime(double time)
{
    std::stringstream ss;
    std::string unit = "ns";
    if(time >= 1000.0)
    {
        time /= 1000.0;
        unit = "us";

        if(time >= 1000.0)
        {
            time /= 1000.0;
            unit = "ms";

            if(time >= 1000.0)
            {
                time /= 1000.0;
                unit = "s";

                if(time >= 60.0)
                {
                    ss << std::floor(time/60.0) << " min " << std::fmod(time, 60.0) << " sec";
                    return ss.str();
                }
            }

        }

    }

    ss << time << ' ' << unit;

    return ss.str();

}

void RandomBenchmark::do_run(size_t nrows, size_t ncols, int num_loops,
                             std::chrono::nanoseconds max_time, bool pretty_print)
{
    std::array<size_t, 2> shape = { nrows, ncols};

    for(auto F : benchmarks)
    {

        xt::random::seed(1);

        std::vector<int64_t> times;
        int64_t times_sum = 0;

        int loop_count;
        for(loop_count = 0; loop_count < num_loops; ++loop_count)
        {
            xt::xtensor<double, 2> elevation = xt::random::rand(shape, 0.0, 1.0);
/*            elevation += 1.0;
            if (nrows > 10)
                for(int i = 3; i<nrows-3; ++i)
                    for(int j = 3; j<ncols-3; ++j)
                        elevation(i,j) -= 1.0;

            for(int i = 0; i<nrows; ++i)
                for(int j = 0; j<ncols; ++j)
                {
                    double x = 1.0 - std::abs((double)j / (double(ncols-1)) *2.0 -1.0);
                    double y = 1.0 - std::abs((double)i / (double(nrows-1)) *2.0 -1.0);

                    elevation(i,j) = std::min(x, y);
                }*/

            xt::xtensor<bool, 2> active_nodes(elevation.shape());
            for(size_t i = 0; i<active_nodes.shape()[0]; ++i)
                for(size_t j = 0; j<active_nodes.shape()[1]; ++j)
                    active_nodes(i,j) = i != 0 && j != 0
                            && i != active_nodes.shape()[0] -1
                            && j != active_nodes.shape()[1] -1;

            auto start = std::chrono::high_resolution_clock::now();

            F.second(elevation, active_nodes);

            auto stop = std::chrono::high_resolution_clock::now();

            times.push_back((stop - start).count());
            times_sum += times.back();

            if (times_sum > max_time.count() && loop_count >= 2)
                break;
        }

        if (loop_count < num_loops)
            ++loop_count;

        double times_avg = (double)times_sum / (double) times.size();
        double time_variance = 0;
        for (int64_t time : times)
            time_variance += ((double)time - times_avg) * ((double)time - times_avg);
        time_variance /= (double) (times.size()-1);

        double time_sdev = std::sqrt(time_variance);

        if (pretty_print)
            std::cout << F.first << ": mean = " << pptime(times_avg)
                      << ", std-dev = " << pptime(time_sdev)
                      << " (" << loop_count << ") loops" << std::endl;
        else
            std::cout << F.first << " " << times_avg
                      << " " << time_sdev << std::endl;
    }
}
