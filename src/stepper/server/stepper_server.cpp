#include <chrono>

#include <hpx/include/iostreams.hpp>
#include <hpx/lcos/gather.hpp>
#include <hpx/lcos/broadcast.hpp>

#include "stepper_server.hpp"


typedef nast_hpx::stepper::server::stepper_server stepper_component;
typedef hpx::components::component<stepper_component> stepper_server_type;

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_COMPONENT(stepper_server_type, stepper_component);
HPX_REGISTER_ACTION(nast_hpx::stepper::server::stepper_server::setup_action,
    stepper_server_setup_action);

typedef std::pair<Real, Real> vec2;
HPX_REGISTER_GATHER(vec2, stepper_server_velocity_gather);

namespace nast_hpx { namespace stepper { namespace server {

stepper_server::stepper_server(uint nl)
: num_localities(nl)
{}

void stepper_server::setup(io::config&& cfg)
{
    // special case for two localities, we want a square configuration
    // of localities
    auto rank = hpx::get_locality_id();

    std::cout << cfg << std::endl;
    std::size_t idx, idy;
    
    
    idx = (rank % cfg.num_localities_x);
    idy = (rank / cfg.num_localities_x);
    part = grid::partition(hpx::find_here(), cfg, idx, idy);
    part.init_sync();
    
    std::size_t step = 0;
    
     std::vector<hpx::naming::id_type> localities;
     for (uint loc = 0; loc < num_localities; loc++)
        localities.push_back(hpx::find_from_basename(stepper_basename, loc).get());
    
    Real dt = cfg.dt;
    for (Real t = 0; ; step++)
    {        
        hpx::future<std::pair<Real, Real> > local_max_velocity =
            part.do_timestep(dt);
    
        // if this is the root locality gather all remote residuals and sum up
        if (hpx::get_locality_id() == 0)
        {
            hpx::future<std::vector<std::pair<Real, Real> > >
            max_velocities =
                hpx::lcos::gather_here(velocity_basename,
                                        std::move(local_max_velocity),
                                        num_localities, step);

            auto local_max_uvs = max_velocities.get();

            std::pair<Real, Real> global_max_uv(0, 0);

            for (auto& max_uv : local_max_uvs)
            {
                global_max_uv.first =
                    (max_uv.first > global_max_uv.first
                        ? max_uv.first : global_max_uv.first);

                global_max_uv.second =
                    (max_uv.second > global_max_uv.second
                        ? max_uv.second : global_max_uv.second);
            }

            Real new_dt =
                std::min(cfg.re / 2. * 1. / (1. / std::pow(cfg.dx, 2)
                            + 1. / std::pow(cfg.dy, 2))
                        ,
                        std::min(cfg.dx / global_max_uv.first,
                                cfg.dy / global_max_uv.second)
                );

            // special case for temperature driven flow
            if (cfg.pr)
                new_dt = std::min(new_dt,
                                (cfg.re * cfg.pr) / 2. * 1.
                                / (1. / std::pow(cfg.dx, 2)
                                + 1. / std::pow(cfg.dy, 2)));

            new_dt *= cfg.tau;

                    

            hpx::lcos::broadcast_apply<set_dt_action>(localities, step, new_dt);
        }
        else
            hpx::lcos::gather_there(velocity_basename, std::move(local_max_velocity),
                                       step).wait();

        if (t >= cfg.t_end)
            break;
            
        dt = dt_buffer.receive(step).get();
        t += dt;
    }
   
}

void stepper_server::set_dt(uint step, Real dt)
{
    dt_buffer.store_received(step, std::move(dt));
}


}//namespace server
}//namespace stepper
}