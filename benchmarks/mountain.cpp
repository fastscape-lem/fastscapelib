#include "dbg_output.hpp"
#include "fastscapelib/basin_graph.hpp"
#include "fastscapelib/flow_routing.hpp"
#include "fastscapelib/bedrock_chanel.hpp"
#include "benchmark_basin_graph.hpp"
#include "benchmark_sinks.hpp"

#include <xtensor/xrandom.hpp>
#include <algorithm>
#include <atomic>
#include <thread>

#include "examples.hpp"
#include <map>

int boruvka_perf = -1;

class FastScapeProgress
{
public:
    FastScapeProgress()
    {
        fast_scape_progress = 0.0;
        print_thread = std::thread(&FastScapeProgress::print_progress, this);
    }
    virtual ~FastScapeProgress()
    {
        fast_scape_progress = 1.0;
        print_thread.join();
    }

    static void update(float p)
    {
        fast_scape_progress = p;
    }

private:
    static std::atomic<double> fast_scape_progress;
    std::thread print_thread;

    void print_progress()
    {
        const int len = 50;

        std::cout << std::unitbuf << "[";
        for (int i = 0; i < len; ++i)
            std::cout << '-';
        std::cout << "]";
        for (int i = 0; i < len + 1; ++i)
            std::cout << '\b';

        int c = 0;
        double progress = 0.0;

        do
        {
            progress = fast_scape_progress.load();

            int nc = int(float(len) * progress);
            while (c < nc)
            {
                std::cout << '#';
                ++c;
            }

            std::this_thread::sleep_for(std::chrono::milliseconds(100));

        } while (progress != 1.0);

        std::cout << std::endl;
    }
};

std::atomic<double> FastScapeProgress::fast_scape_progress(0.0);

template <class Elev_T, class Active_T, class Uplift_XT>
std::vector<std::pair<index_t, double>>
fastscape(Elev_T& elevation,
          const Active_T& active_nodes,
          double dx,
          double dt,
          int num_iter,
          Uplift_XT&& uplift,
          bool correct_flow = true)
{
    auto shape = elevation.shape();
    int nrows = (int) shape[0];
    int ncols = (int) shape[1];
    std::array<size_t, 1> shape1D = { (size_t) nrows * (size_t) ncols };

    double dy = dx;

    fastscapelib::BasinGraph<index_t, index_t, double> bg;
    xt::xtensor<index_t, 1> basins(shape1D);
    xt::xtensor<index_t, 1> stack(shape1D);
    xt::xtensor<index_t, 1> receivers(shape1D);
    xt::xtensor<double, 1> dist2receivers(shape1D);
    xt::xtensor<index_t, 1> ndonors(shape1D);
    xt::xtensor<index_t, 2> donors({ shape1D[0], 8 });
    xt::xtensor<double, 1> area(shape1D);
    xt::xtensor<double, 1> erosion(shape1D);

    std::vector<std::pair<index_t, double>> out;


    for (int s = 0; s < num_iter; ++s)
    {
        FastScapeProgress::update(float(s) / float(num_iter - 1));

        // uplift
        elevation += xt::cast<double>(active_nodes) * uplift * dt;

        fastscapelib::compute_receivers_d8(
            receivers, dist2receivers, elevation, active_nodes, dx, dy);

        fastscapelib::compute_donors(ndonors, donors, receivers);
        fastscapelib::compute_stack(stack, ndonors, donors, receivers);

        if (correct_flow)
            fastscapelib::correct_flowrouting<fs::BasinAlgo::Boruvka, fs::ConnectType::Carved>(
                bg,
                basins,
                receivers,
                dist2receivers,
                ndonors,
                donors,
                stack,
                active_nodes,
                elevation,
                dx,
                dy);

        fs::compute_drainage_area(area, stack, receivers, dx, dy);
        fs::erode_spower(erosion,
                         elevation,
                         stack,
                         receivers,
                         dist2receivers,
                         area,
                         7.0e-4,
                         0.4,
                         1.0,
                         dt,
                         1.0e-4);

        for (size_t k = 0; k < nrows * ncols; ++k)
            elevation(k) = elevation(k) - erosion(k);

        out.emplace_back(bg.basin_count(), xt::mean(elevation)[0]);
    }

    return out;
}

template <class T>
typename std::decay_t<T>::value_type
fetch(const T& a, double r, double c)
{
    using F = typename std::decay_t<T>::value_type;

    auto shape = a.shape();

    // r \in [0, nr-1]

    int ir = std::min((int) r, (int) shape[0] - 2);
    int ic = std::min((int) c, (int) shape[1] - 2);

    double sr = r - double(ir);
    double sc = c - double(ic);

    F a00 = a(ir, ic);
    F a01 = a(ir, ic + 1);
    F a10 = a(ir + 1, ic);
    F a11 = a(ir + 1, ic + 1);

    return a00 * (1.0 - sr) * (1.0 - sc) + a01 * (1.0 - sr) * sc + a10 * sr * (1.0 - sc)
           + a11 * sc * sr;
}

template <class T>
void
sample(const T& in, T& out)
{
    for (int r = 0; r < out.shape()[0]; ++r)
        for (int c = 0; c < out.shape()[1]; ++c)
        {
            double dr = double(r) / double(out.shape()[0] - 1) * double(in.shape()[0] - 1);
            double dc = double(c) / double(out.shape()[1] - 1) * double(in.shape()[1] - 1);
            out(r, c) = fetch(in, dr, dc);
        }
}

void
generate_mountain()
{
    int n = 32;
    double dx = 64 * 10;
    double dt = 1000000;


    xt::xtensor<double, 2> elevation = xt::random::rand({ (size_t) n, (size_t) n }, 0.0, 1e-3);

    while (n <= 4096 * 4)
    // while (n <= 0)
    {
        xt::xtensor<bool, 2> active_nodes(elevation.shape());

        for (size_t i = 0; i < active_nodes.shape()[0]; ++i)
            for (size_t j = 0; j < active_nodes.shape()[1]; ++j)
                active_nodes(i, j) = i != 0 && j != 0 && i != active_nodes.shape()[0] - 1
                                     && j != active_nodes.shape()[1] - 1;

        int num_iter = n == 32 ? 100 : 4;

        fastscape(elevation, active_nodes, dx, dt, num_iter, 5e-3);

        dbg_out("results/mountain/bedrock-", n, elevation, elevation.shape());

        xt::xtensor<double, 2> tmp({ (size_t) (2 * n), (size_t) (2 * n) }, -1);

        // linear upsampling
        sample(elevation, tmp);

        elevation = tmp;
        n *= 2;
        dx *= 0.5;
        dt *= 0.5;
    }
}

template <class S>
std::vector<size_t>
holes(const S& shape)
{
    std::vector<size_t> ids;
    ids.reserve((shape[0] - 1) / 2 * ((shape[1] - 1) / 2));
    for (size_t r = 1; r < shape[0] - 1; r += 2)
        for (size_t c = 1; c < shape[1] - 1; c += 2)
            ids.push_back(c + r * shape[1]);

    assert(ids.size() == (shape[0] - 1) / 2 * ((shape[1] - 1) / 2));

    std::mt19937 gen = std::mt19937(123);  // seed = 123
    std::shuffle(ids.begin(), ids.end(), gen);

    return ids;
}

template <class T>
void
dig_hole(T& a, size_t hid)
{
    size_t r = hid / a.shape()[1];
    size_t c = hid % a.shape()[1];

    double h = std::numeric_limits<double>::max();
    for (size_t rr = r - 1; rr <= r + 1; ++rr)
        for (size_t cc = c - 1; cc <= c + 1; ++cc)
            h = std::min(h, a(rr, cc));

    a(r, c) = h - 1e-5;
}


void
example_mountain()
{
    double prop_multiplier = 0.1;

    xt::xtensor<double, 1> h_prop = { 0,       0.0078125,
                                      .015625, 0.5 * (.015625 + .03125),
                                      .03125,  0.5 * (.03125 + .0625),
                                      .0625,   0.5 * (.0625 + .125),
                                      .125,    0.5 * (.125 + .25),
                                      .25,     0.5 * (.25 + .5),
                                      .5,      0.5 * (.5 + 1.0),
                                      1.0 };
    xt::xtensor<int, 1> m_size
        = { /*32, 64, 128, 256, 512, 1024, 2048,*/ 4096, 4096 * 2, 4096 * 4 };
    // xt::xtensor<int, 1> m_size = { 32, 64, 128, 256 };

    xt::xtensor<double, 2> results({ m_size.size(), h_prop.size() });

    xt::xtensor<double, 2> elevation;

    using FastscapeFunctionType = std::function<void(xt::xtensor<bm_index, 1>& /*stack*/,
                                                     xt::xtensor<bm_index, 1>& /*receivers*/,
                                                     xt::xtensor<double, 1>& /*dist2receviers*/,
                                                     const xt::xtensor<double, 2>& /*elevation*/,
                                                     const xt::xtensor<bool, 2>& /*active nodes  */,
                                                     double dx,
                                                     double dy)>;

    std::map<std::string, FastscapeFunctionType> funcs;
    /*funcs["Kruskal simple"] = benchmark_fastscape_basin<fs::BasinAlgo::Kruskal,
    fs::ConnectType::Simple>; funcs["Boruvka simple"] =
    benchmark_fastscape_basin<fs::BasinAlgo::Boruvka, fs::ConnectType::Simple>; funcs["Kruskal
    carved"] = benchmark_fastscape_basin<fs::BasinAlgo::Kruskal, fs::ConnectType::Carved>;
    funcs["Boruvka carved"] = benchmark_fastscape_basin<fs::BasinAlgo::Boruvka,
    fs::ConnectType::Carved>;*/
    funcs["Kruskal sloped"]
        = benchmark_fastscape_basin<fs::BasinAlgo::Kruskal, fs::ConnectType::Sloped>;
    funcs["Boruvka sloped"]
        = benchmark_fastscape_basin<fs::BasinAlgo::Boruvka, fs::ConnectType::Sloped>;
    funcs["Sinks"] = benchmark_fastscape_sinks;
#ifdef ENABLE_RICHDEM
    funcs["Wei2018"] = benchmark_fastscape_wei2018;
    funcs["Zhou2016"] = benchmark_fastscape_zhou2016;
#endif


    for (int s : m_size)
    {
        std::stringstream out;
        out << "bench = {";

        out << s << ":{";

        dbg_in("results/mountain/bedrock-", s, elevation);

        auto holes_id = holes(elevation.shape());

        int num_holes = 0;
        for (double hp : h_prop)
        {
            int target_holes = (int) (prop_multiplier * hp * double(holes_id.size()));
            int max_holes
                = (int) (prop_multiplier * h_prop[h_prop.size() - 1] * double(holes_id.size()));

            out << hp << ":{'numholes':" << target_holes << ',';

            for (int i = num_holes; i < target_holes; ++i)
                dig_hole(elevation, holes_id[i]);
            num_holes = target_holes;


            std::array<size_t, 1> shape1D = { elevation.size() };

            xt::xtensor<bm_index, 1> stack(shape1D);
            xt::xtensor<bm_index, 1> receivers(shape1D);
            xt::xtensor<double, 1> dist2receivers(shape1D);
            xt::xtensor<bool, 2> active_nodes(elevation.shape());
            for (size_t i = 0; i < active_nodes.shape()[0]; ++i)
                for (size_t j = 0; j < active_nodes.shape()[1]; ++j)
                    active_nodes(i, j) = i != 0 && j != 0 && i != active_nodes.shape()[0] - 1
                                         && j != active_nodes.shape()[1] - 1;

            for (auto f : funcs)
            {
                std::vector<double> times;
                double times_sum = 0;

                boruvka_perf = -1;

                for (int i = 0; i < 10; ++i)
                {
                    std::cerr << s << ' ' << num_holes << '/' << max_holes << ' ' << f.first << ' '
                              << i + 1 << "/10" << std::endl;

                    xt::xtensor<double, 2> elev = elevation;

                    auto start = std::chrono::high_resolution_clock::now();

                    f.second(stack, receivers, dist2receivers, elev, active_nodes, 10.0, 10.0);

                    auto stop = std::chrono::high_resolution_clock::now();

                    std::chrono::duration<double, std::milli> fp_ms = stop - start;

                    times.push_back(fp_ms.count());
                    times_sum += times.back();
                }
                double times_avg = times_sum / (double) times.size();
                double time_variance = 0;
                for (double time : times)
                    time_variance += (time - times_avg) * (time - times_avg);
                time_variance /= (double) (times.size() - 1);

                double time_sdev = std::sqrt(time_variance);

                out << "'" << f.first << "':(" << times_avg << ',' << time_sdev << ","
                    << boruvka_perf << "),";
            }

            out << "},";  // hp
        }
        out << "},";  // s

        out << "}";

        std::stringstream filename;
        filename << "results/mountain/bench" << s << ".py";
        std::ofstream file(filename.str());
        if (!file)
            std::cerr << "Impossible to open file " << filename.str() << std::endl;
        else
            file << out.str();
    }
}

template <class Shape>
xt::xtensor<double, 2>
simple_noise(Shape shape, int frequency)
{
    xt::xtensor<double, 2> out(shape),
        noise = xt::random::rand({ frequency + 1, frequency + 1 }, 0.0, 1.0);

    sample(noise, out);
    return out;
}


void
fastscape_pits()
{
    for (int k = 0; k < 2; ++k)
        for (int n = 1024 * 4; n <= 1024 * 16; n *= 2)
        {
            double dx = 100;
            double dt = 10000;
            int niter = 20;

            if (k == 1)
                dx *= 1024 * 4 / double(n);

            xt::xtensor<double, 1> noise_ampl = xt::linspace(1e-10, 0.1, 4);

            xt::xtensor<double, 2> elevation_init
                = xt::random::rand({ (size_t) n, (size_t) n }, 0.0, 1e-1 * dt);

            xt::xtensor<bool, 2> active_nodes(elevation_init.shape());

            for (size_t i = 0; i < active_nodes.shape()[0]; ++i)
                for (size_t j = 0; j < active_nodes.shape()[1]; ++j)
                    active_nodes(i, j) = i != 0 && j != 0 && i != active_nodes.shape()[0] - 1
                                         && j != active_nodes.shape()[1] - 1;

            std::stringstream out;
            out << "pits = {";

            size_t c = 0;
            for (auto ampl : noise_ampl)
            {
                // xt::xtensor<double, 2> elevation = elevation_init;

                xt::xtensor<double, 2> elevation = xt::zeros_like(active_nodes);

                std::cout << "Pits: " << 100 * c++ / (noise_ampl.size() - 1) << "% ";
                FastScapeProgress progress;

                xt::xtensor<double, 2> uplift
                    = ((2.0 * simple_noise(elevation.shape(), 1500) - 1.0) * ampl + 1.0) * 5e-3;
                auto basin_count = fastscape(elevation, active_nodes, dx, dt, niter, uplift);

                out << ampl << ": [";
                std::for_each(basin_count.begin(),
                              basin_count.end(),
                              [&out](auto c) { out << c.first << ","; });
                out << "],";
            }
            out << "}";

            std::stringstream filename;
            filename << "results/mountain/pits" << n << '-' << k << ".py";
            std::ofstream file(filename.str());
            if (!file)
                std::cerr << "Impossible to open file results/mountain/pits.py\n";
            else
                file << out.str();
        }
}

void
escarpment()
{
    int nx = 1024;
    int ny = 256;
    double dx = 100;
    double dt = 1000;
    int niter = 70;

    double global_height = 500;


    xt::xtensor<double, 2> elevation
        = xt::random::rand({ (size_t) ny, (size_t) nx }, 0.0, 1e-3) + global_height;

    xt::xtensor<bool, 2> active_nodes(elevation.shape(), true);

    for (size_t i = 0; i < active_nodes.shape()[0]; ++i)
    {
        active_nodes(i, 0) = false;
        elevation(i, 0) = 0.0;
    }


    xt::xtensor<double, 2> elev0 = elevation;
    {
        std::cout << "Loc min: false ";
        FastScapeProgress progress;
        fastscape(elevation, active_nodes, dx, dt, niter, 0.0, false);
    }
    {
        std::cout << "Loc min: true ";
        FastScapeProgress progress;
        fastscape(elev0, active_nodes, dx, dt, niter, 0.0);
    }

    dbg_out("results/mountain/escarp", 0, elevation, elevation.shape());
    dbg_out("results/mountain/escarp", 1, elev0, elev0.shape());
}


void
locmin2()
{
    int nx = 500;
    int ny = 500;
    double dx = 100;
    double dt = 1000;
    int niter = 500;

    double global_height = 0.0;


    xt::xtensor<double, 2> elevation
        = xt::random::rand({ (size_t) ny, (size_t) nx }, 0.0, 1e-3) + global_height;

    xt::xtensor<bool, 2> active_nodes(elevation.shape(), true);

    for (size_t i = 0; i < active_nodes.shape()[0]; ++i)
        for (size_t j = 0; j < active_nodes.shape()[1]; ++j)
            active_nodes(i, j) = i != 0 && j != 0 && i != active_nodes.shape()[0] - 1
                                 && j != active_nodes.shape()[1] - 1;

    std::stringstream out;
    out << "ml = {False:[";

    std::stringstream out2;
    out2 << "mlp = {False:[";


    xt::xtensor<double, 2> elev0 = elevation;
    {
        std::cout << "Loc min: false ";
        FastScapeProgress progress;
        auto mean_height = fastscape(elevation, active_nodes, dx, dt, niter, 5e-3, false);

        std::for_each(
            mean_height.begin(), mean_height.end(), [&out](auto c) { out << c.second << ","; });
        std::for_each(
            mean_height.begin(), mean_height.end(), [&out2](auto c) { out2 << c.first << ","; });
    }
    out << "],True:[";
    out2 << "],True:[";
    {
        std::cout << "Loc min: true ";
        FastScapeProgress progress;
        auto mean_height = fastscape(elev0, active_nodes, dx, dt, niter, 5e-3);

        std::for_each(
            mean_height.begin(), mean_height.end(), [&out](auto c) { out << c.second << ","; });
        std::for_each(
            mean_height.begin(), mean_height.end(), [&out2](auto c) { out2 << c.first << ","; });
    }
    out << "]}";
    out2 << "]}";

    std::ofstream file("results/mountain/ml.py");
    if (!file)
        std::cerr << "Impossible to open file results/mountain/ml.py\n";
    else
        file << out.str() << "\n" << out2.str();

    // dbg_out("results/mountain/escarp", 0, elevation, elevation.shape());
    // dbg_out("results/mountain/escarp", 1, elev0, elev0.shape());
}
