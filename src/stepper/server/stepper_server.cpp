#include <chrono>

//#include <hpx/include/iostreams.hpp>
#include <hpx/lcos/gather.hpp>
#include <hpx/lcos/broadcast.hpp>
#include <hpx/runtime/components/migrate_component.hpp>

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
    auto rank = hpx::get_locality_id();
    auto num_localities = hpx::get_num_localities_sync();

    cfg.idx = (rank % cfg.num_localities_x);
    cfg.idy = (rank / cfg.num_localities_x);

    cfg.num_partitions_x = cfg.num_localities_x;
    cfg.num_partitions_y = cfg.num_localities_y;
    cfg.num_partitions = cfg.num_localities;

    Real dt = cfg.initial_dt;
    Real dx = cfg.dx;
    Real dy = cfg.dy;
    Real re = cfg.re;
    Real pr = cfg.pr;
    Real tau = cfg.tau;
    Real t_end = cfg.t_end;

    std::size_t max_timesteps = cfg.max_timesteps;

    grid::partition part(hpx::find_here(), std::move(cfg));
    part.init_sync();

    std::vector<hpx::naming::id_type> localities;
    for (uint loc = 0; loc < num_localities; loc++)
        localities.push_back(hpx::find_from_basename(stepper_basename, loc).get());

    std::size_t step = 0;
    bool running = true;

    for (Real t = 0;; step++)
    {
        if (max_timesteps > 0 && step >= max_timesteps)
            break;
       // std::cout << "step " << step << std::endl;
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

            max_velocities.then(
                hpx::util::unwrapped(
                    [=, &t](std::vector<std::pair<Real, Real> > local_max_uvs)
                    {
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
                            std::min(re / 2. * 1. / (1. / std::pow(dx, 2)
                                        + 1. / std::pow(dy, 2))
                                    ,
                                    std::min(dx / global_max_uv.first,
                                            dy / global_max_uv.second)
                            );

                        new_dt *= tau;

                        hpx::lcos::broadcast_apply<set_dt_action>(localities, step, new_dt);
                    }
                )
            );
        }
        else
            hpx::lcos::gather_there(velocity_basename, std::move(local_max_velocity),
                                       step);

        //std::cout << "receiving in step " << step << std::endl;
        if (t >= t_end)
            break;
        t += dt;
        dt = dt_buffer.receive(step).get();
    }
}

void stepper_server::set_dt(uint step, Real dt)
{
    dt_buffer.store_received(step, std::move(dt));
}


}//namespace server
}//namespace stepper
}
