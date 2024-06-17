#include <algorithm>

#include "xtensor/xtensor.hpp"
#include "xtensor/xstrided_view.hpp"

#include "gtest/gtest.h"

#include "fastscapelib/flow/flow_graph.hpp"
#include "fastscapelib/flow/flow_router.hpp"
#include "fastscapelib/flow/flow_snapshot.hpp"
#include "fastscapelib/grid/raster_grid.hpp"


namespace fs = fastscapelib;

namespace fastscapelib
{
    namespace testing
    {
        template <class T>
        struct Scalar
        {
            using type = T;
        };

        template <class T>
        struct Vector
        {
            using type = xt::xtensor<T, 1>;
        };

        template <class T>
        struct Grid
        {
            using type = xt::xtensor<T, 2>;
        };

        template <class T>
        using scalar_t = typename Scalar<T>::type;

        template <class T>
        using vector_t = typename Vector<T>::type;

        template <class T>
        using grid_t = typename Grid<T>::type;

        template <template <typename> class Selector>
        struct Data
        {
            Selector<double> erosion;
            Selector<double> drainage_area;
        };

        struct ReceiversData : public Data<vector_t>
        {
            vector_t<double> weight;
            vector_t<double> distance;
        };

        struct SplEroderKernelNodeData : public Data<scalar_t>
        {
            ReceiversData receivers;
        };

        using SplEroderKernelData = Data<grid_t>;

        class SplEroder
        {
        private:
            static int kernel_func(void* node_data_ptr)
            {
                auto& node = *reinterpret_cast<SplEroderKernelNodeData*>(node_data_ptr);

                return 0;
            }

            static int ndata_getter(std::size_t idx, void* data_ptr, void* node_data_ptr)
            {
                auto& data = *reinterpret_cast<SplEroderKernelData*>(data_ptr);
                auto& node_data = *reinterpret_cast<SplEroderKernelNodeData*>(node_data_ptr);

                node_data.erosion = data.erosion[idx];

                return 0;
            }

            static int ndata_setter(std::size_t idx, void* node_data_ptr, void* data_ptr)
            {
                auto& node_data = *reinterpret_cast<SplEroderKernelNodeData*>(node_data_ptr);
                auto& data = *reinterpret_cast<SplEroderKernelData*>(data_ptr);

                data.erosion[idx] = node_data.erosion;

                return 0;
            }

            static void* ndata_create()
            {
                SplEroderKernelNodeData* node_data_ptr = new SplEroderKernelNodeData;

                return reinterpret_cast<void*>(node_data_ptr);
            }

        public:
            using array_type = typename xt::xtensor<double, 2>;

            template <class FG>
            SplEroder(FG& /*flow_graph*/, int threads_count = 1)
            {
                kernel.func = &SplEroder::kernel_func;
                kernel.node_data_getter = &SplEroder::ndata_getter;
                kernel.node_data_setter = &SplEroder::ndata_setter;
                kernel.node_data_create = &SplEroder::ndata_create;
                kernel.node_data_init = nullptr;
                kernel.node_data_free = [](void* node_data_ptr)
                { delete reinterpret_cast<SplEroderKernelNodeData*>(node_data_ptr); };
                kernel.n_threads = threads_count;
                kernel.apply_dir = flow_graph_traversal_dir::breadth_upstream;

                kernel_data.data = reinterpret_cast<void*>(&data);
            }

            array_type erode(array_type& elevation, double /*dt*/)
            {
                return elevation;
            }

            flow_kernel kernel;
            flow_kernel_data kernel_data;
            SplEroderKernelData data;
        };

        TEST(flow_kernel, spl_eroder)
        {
            using flow_graph_type = fs::flow_graph<fs::raster_grid<>>;
            using grid_type = fs::raster_grid<>;

            grid_type grid = grid_type({ 4, 4 }, { 1.0, 1.0 }, fs::node_status::fixed_value);

            auto graph = flow_graph_type(grid, { fs::single_flow_router() });

            xt::xtensor<double, 2> elevation{ { 0.6, 0.6, 0.6, 0.6 },
                                              { 0.4, 0.4, 0.4, 0.4 },
                                              { 0.2, 0.2, 0.2, 0.2 },
                                              { 0.1, 0.0, 0.1, 0.1 } };

            SplEroder eroder(graph);
            eroder.data.erosion = xt::zeros_like(elevation);

            graph.apply_kernel(eroder.kernel, eroder.kernel_data);
            eroder.erode(elevation, 2e4);
        }
    }
}
