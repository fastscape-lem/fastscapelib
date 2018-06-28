#include "random_benchmark.hpp"
#include <chrono>

#include <xtensor/xtensor.hpp>
#include <xtensor/xrandom.hpp>

#include "examples.hpp"


using namespace std::literals;


int main(int , char* [])
{
	//example_vornoi();
	//example_jail();
	//generate_mountain();
	example_mountain();

	std::system("pause");

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
