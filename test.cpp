#include <iostream>
#include <cassert>

#include "xtensor/xio.hpp"

#include "include/fastscape/fastscape.hpp"


namespace fs = fastscape;


void test() {
    xt::xtensor<double, 2> arr1
       {{1.0, 2.0, 3.0},
        {2.0, 0.1, 7.0},
        {2.0, 5.0, 7.0}};

    std::cout << arr1 << std::endl;

    fs::fill_sinks_flat(arr1);

    std::cout << arr1 << std::endl;

    xt::xtensor<double, 2> arr2
       {{1.00001, 2.0, 3.0},
        {2.0,     0.1, 7.0},
        {2.0,     5.0, 7.0}};

    std::cout << arr2 << std::endl;

    assert(arr1(1, 1) == arr1(0, 0));

    fs::fill_sinks_sloped(arr2);

    std::cout << arr2 << std::endl;

    assert(arr2(1, 1) > arr2(0, 0));

    std::cout << "end" << std::endl;
}
