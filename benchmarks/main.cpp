#include "random_benchmark.hpp"
#include <chrono>

#include <xtensor/xtensor.hpp>
#include <xtensor/xrandom.hpp>

#include "examples.hpp"


using namespace std::literals;


int
main(int i, char* p[])
{
    /*
    auto now = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 100000; ++i)
        volatile int a = i;
    auto tmp = std::chrono::high_resolution_clock::now();

    double ticks = (tmp - now).count();
    double secs = std::chrono::duration_cast<std::chrono::duration<double>>(tmp - now).count();

    std::cout << 10.0e10 / ticks * secs << std::endl;


        std::cin.get();
    return 0;*/

    std::string param = i > 1 ? p[1] : "fastscape_pits";

    if (param == "vornoi" || param == "all")
        example_vornoi();
    if (param == "jail" || param == "all")
        example_jail();
    if (param == "gen_mountain" || param == "all")
        generate_mountain();
    if (param == "mountain" || param == "all")
        example_mountain();
    if (param == "fastscape_pits" || param == "all")
        fastscape_pits();
    if (param == "escarpment" || param == "all")
        escarpment();
    if (param == "locmin2" || param == "all")
        locmin2();

    std::cout << "Press enter to continue ...";
    std::cin.get();

    return 0;

    // for (size_t x = 1; x<50000; x*=2)
    // for (size_t x = 1; x<5000; x+=1)
    // for (size_t x = 4096; x<4097; x*=2)
    // for (size_t x = 16*1024; x<=16*1024; x*=2)
    // for (size_t x = 64; x <= 8 * 1024; x +=64) //128 iter max
    for (size_t x = 64; x <= 4097; x *= 2)
    {
        std::cout << "Size " << x << std::endl;
        Benchmark::run("fastscape", x, x, 3, 2s, true);
    }
    char c;
    std::cin >> c;
}
