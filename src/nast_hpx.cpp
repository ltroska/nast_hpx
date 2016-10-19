#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include "io/config.hpp"
#include "grid/partition.hpp"

int hpx_main(boost::program_options::variables_map& vm)
{
    const auto cfg_path = vm["cfg"].as<std::string>();
    const auto iterations = vm["iterations"].as<std::size_t>();
    const auto timesteps = vm["timesteps"].as<std::size_t>();

    nast_hpx::io::config cfg = nast_hpx::io::config::read_config_from_file(cfg_path.c_str(), hpx::get_locality_id(), 1);
    cfg.max_timesteps = timesteps;
    cfg.verbose = vm.count("verbose") ? true : false;

    std::cout << cfg << std::endl;

    if (cfg.verbose)
        std::cout << "Threads on locality " << hpx::get_locality_id()
            << " = " << hpx::get_os_thread_count() << std::endl;

    std::cout
        << "Running simulation on " << cfg.i_max + 2 << "x" << cfg.j_max + 2
        << " cells on " << 1 << " nodes ";
        if (timesteps == 0)
            std::cout << "until t_end " << cfg.t_end << std::endl;
        else
            std::cout << "for " << timesteps << " iterations"
        << " and " << iterations << std::endl;

    nast_hpx::grid::partition part(hpx::find_here(), cfg);
    part.init_sync();

    Real t = 0;
    Real dt = cfg.initial_dt;

    std::size_t step = 0;

    while (t < cfg.t_end)
    {
        if (cfg.max_timesteps > 0 && step >= cfg.max_timesteps)
            break;

        hpx::future<nast_hpx::pair<Real> > max_velocity =
         part.do_timestep(dt);

        t += dt;

        dt = max_velocity.then(
          hpx::util::unwrapped(
              [&cfg](nast_hpx::pair<Real> max_uv)
              -> Real
              {
                  Real new_dt =
                      std::min(cfg.re / 2. * 1. / (1. / std::pow(cfg.dx, 2)
                                  + 1. / std::pow(cfg.dy, 2))
                              ,
                              std::min(cfg.dx / max_uv.x,
                                      cfg.dy / max_uv.y)
                      );

                  return new_dt * cfg.tau;
              }
        )).get();

      ++step;
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
