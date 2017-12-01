#include <iostream>

#include "xtensor/xio.hpp"

#include "depressions.hpp"


void test() {
    xt::xtensor<double, 2> arr1
       {{1.0, 2.0, 3.0},
        {2.0, 0.1, 7.0},
        {2.0, 5.0, 7.0}};

    std::cout << arr1 << std::endl;

    priority_flood_original(arr1);

    std::cout << arr1 << std::endl;
}
