#include "random_benchmark.hpp"

#include "xtensor/xrandom.hpp"

void random_run(size_t rows, size_t cols, RandomFunctionType func)
{
    std::array<size_t, 2> shape = { rows, cols};

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

    func(elevation, active_nodes);
}
