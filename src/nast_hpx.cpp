#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include "io/config.hpp"
#include "stepper/stepper.hpp"

int hpx_main(boost::program_options::variables_map& vm)
{
    const auto cfg_path = vm["cfg"].as<std::string>();
    const auto iterations = vm["iterations"].as<std::size_t>();
    const auto timesteps = vm["timesteps"].as<std::size_t>();

    nast_hpx::io::config cfg = nast_hpx::io::config::read_config_from_file(cfg_path.c_str(), hpx::get_locality_id(), hpx::get_num_localities_sync());
    cfg.max_timesteps = timesteps;
    cfg.verbose = vm.count("verbose") ? true : false;

    if (cfg.verbose)
        std::cout << "Threads on locality " << hpx::get_locality_id()
            << " = " << hpx::get_os_thread_count() << std::endl;

    std::cout
        << "Running simulation on " << cfg.i_max + 2 << "x" << cfg.j_max + 2
        << " cells on " << cfg.num_localities << " nodes ";

    if (timesteps == 0)
        std::cout << "until t_end " << cfg.t_end << std::endl;
    else
        std::cout << "for " << timesteps << " iterations"
            << " and " << iterations << " runs!" << std::endl;

    auto rank = hpx::get_locality_id();
    auto num_localities = hpx::get_num_localities_sync();

    cfg.idx = (rank % (cfg.num_localities_x * cfg.num_localities_y)) % cfg.num_localities_x;
    cfg.idy = (rank % (cfg.num_localities_x * cfg.num_localities_y)) / cfg.num_localities_x;
    cfg.idz = rank / (cfg.num_localities_x * cfg.num_localities_y);

    cfg.num_partitions_x = cfg.num_localities_x;
    cfg.num_partitions_y = cfg.num_localities_y;
    cfg.num_partitions_z = cfg.num_localities_z;
    cfg.num_partitions = cfg.num_localities;

    double avgtime = 0.;
    double maxtime = 0.;
    double mintime = 365. * 24. * 3600.;

    std::vector<hpx::performance_counters::performance_counter> idle_rate_counters(num_localities);
    std::vector<std::size_t> avg_idle_rates(num_localities, 0);
    std::vector<std::size_t> max_idle_rates(num_localities, 0);
    std::vector<std::size_t> min_idle_rates(num_localities, 20000);
    std::size_t avgidlerate = 0;
    std::size_t maxidlerate = 0;
    std::size_t minidlerate = 20000;

    for (std::size_t loc = 0; loc < num_localities; ++loc)
        idle_rate_counters[loc] = hpx::performance_counters::performance_counter("/threads{locality#" + std::to_string(loc) + "/total}/idle-rate");

    nast_hpx::stepper::stepper step;
    step.setup(cfg);

    for (std::size_t iter = 0; iter < iterations; ++iter)
    {
        hpx::util::high_resolution_timer t;

        step.run();

        double elapsed = t.elapsed();

        if (iter > 0 || iterations == 1)
        {
            avgtime += elapsed;
            maxtime = std::max(maxtime, elapsed);
            mintime = std::min(mintime, elapsed);

            for (std::size_t loc = 0; loc < num_localities; ++loc)
            {
                std::size_t idle_rate = idle_rate_counters[loc].get_value<std::size_t>().get();

                avg_idle_rates[loc] += idle_rate;
                max_idle_rates[loc] = std::max(max_idle_rates[loc], idle_rate);
                min_idle_rates[loc] = std::min(min_idle_rates[loc], idle_rate);

                avgidlerate += idle_rate;
                maxidlerate = std::max(maxidlerate, idle_rate);
                minidlerate = std::min(minidlerate, idle_rate);
            }
        }
    }

    if (rank == 0)
    {
        avgtime = avgtime / static_cast<double>(
                    (std::max)(iterations-1, static_cast<boost::uint64_t>(1)));

        avgidlerate = avgidlerate / static_cast<double>(
                    (std::max)(iterations-1, static_cast<boost::uint64_t>(1)))
                    / num_localities;

        for (std::size_t loc = 0; loc < num_localities; ++loc)
            avg_idle_rates[loc] = avg_idle_rates[loc] / static_cast<double>(
                    (std::max)(iterations-1, static_cast<boost::uint64_t>(1)));

        std::cout
            << "Avg time (s):\t" << avgtime << "\n"
            << "Min time (s):\t" << mintime << "\n"
            << "Max time (s):\t" << maxtime << "\n"
            << "Avg idle rate:\t" << avgidlerate << "\n"
            << "Min idle rate:\t" << minidlerate << "\n"
            << "Max idle rate:\t" << maxidlerate << "\n";

        for (std::size_t loc = 0; loc < num_localities; ++loc)
           std::cout
            << "Avg idle rate (" << loc << "):\t" << avg_idle_rates[loc] << "\n"
            << "Min idle rate (" << loc << "):\t" << min_idle_rates[loc] << "\n"
            << "Max idle rate (" << loc << "):\t" << max_idle_rates[loc] << "\n";

        std::cout << std::endl;

    }


    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    using namespace boost::program_options;

    options_description desc_commandline;

    desc_commandline.add_options()
    ("cfg", value<std::string>()->required(),
         "path to config xml")
    ("iterations", value<std::size_t>()->default_value(1),
         "Number of runs of the simulation")
    ("timesteps", value<std::size_t>()->default_value(0),
         "Number of timesteps per run (0 = use t_end)")
     ( "verbose", "Verbose output");

    std::vector<std::string> cfg;
    cfg.push_back("hpx.run_hpx_main!=1");

    return hpx::init(desc_commandline, argc, argv, cfg);
}
