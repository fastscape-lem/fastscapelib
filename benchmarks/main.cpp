#include "random_benchmark.hpp"
#include <chrono>

#include <xtensor/xtensor.hpp>
#include <xtensor/xrandom.hpp>

#include "examples.hpp"


using namespace std::literals;


int main(int i, char* p[])
{

	std::string param = i > 1 ? p[1] : "escarpment";

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

	std::cout << "Press enter to continue ...";
	std::cin.get();

	return 0;

    //for (size_t x = 1; x<50000; x*=2)
    //for (size_t x = 1; x<5000; x+=1)
    //for (size_t x = 4096; x<4097; x*=2)
    //for (size_t x = 16*1024; x<=16*1024; x*=2)
    //for (size_t x = 64; x <= 8 * 1024; x +=64) //128 iter max
    for (size_t x = 64; x<= 4097; x*=2)
    {
        std::cout << "Size " << x << std::endl;
        Benchmark::run("fastscape", x, x, 3, 2s, true);
    }
	char c;
	std::cin >> c;
}
