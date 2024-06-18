#include "fastscapelib/grid/healpix_grid.hpp"

#include "gtest/gtest.h"


namespace fs = fastscapelib;


namespace fastscapelib
{
    namespace testing
    {

        TEST(healpix_grid, ctor)
        {
            auto grid = fs::healpix_grid<>(25);
        }
    }
}
