#include "random_benchmark.hpp"

int main(int , char* [])
{

    for (size_t x = 1; x<5000; x*=2)
    {
        std::cout << "Size " << x << std::endl;
        RandomBenchmark::run(x, x, 200, true);
    }

}
