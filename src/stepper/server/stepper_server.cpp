#include <chrono>

//#include <hpx/include/iostreams.hpp>
#include <hpx/lcos/gather.hpp>
#include <hpx/lcos/broadcast.hpp>
#include <hpx/runtime/components/migrate_component.hpp>

#include "stepper_server.hpp"

#include "util/triple.hpp"


typedef nast_hpx::stepper::server::stepper_server stepper_component;
typedef hpx::components::component<stepper_component> stepper_server_type;

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_COMPONENT(stepper_server_type, stepper_component);
HPX_REGISTER_ACTION(nast_hpx::stepper::server::stepper_server::setup_action,
    stepper_server_setup_action);
HPX_REGISTER_ACTION(nast_hpx::stepper::server::stepper_server::run_action,
    stepper_server_run_action);

typedef nast_hpx::triple<Real> vec3;
HPX_REGISTER_GATHER(vec3, stepper_server_velocity_gather);

namespace nast_hpx { namespace stepper { namespace server {

stepper_server::stepper_server(uint nl)
: num_localities(nl)
{}

void stepper_server::setup(io::config const& cfg)
{
    rank = hpx::get_locality_id();
    num_localities = hpx::get_num_localities_sync();

    dx = cfg.dx;
    dy = cfg.dy;
    dz = cfg.dz;
    re = cfg.re;
    pr = cfg.pr;
    tau = cfg.tau;
    t_end = cfg.t_end;
    init_dt = cfg.initial_dt;



    max_timesteps = cfg.max_timesteps;
    step = 0;

    part = grid::partition(hpx::find_here(), cfg);

    for (uint loc = 0; loc < num_localities; loc++)
        localities.push_back(hpx::find_from_basename(stepper_basename, loc).get());
}

void stepper_server::run()
{
    part.init_sync();

    Real dt = init_dt;

    std::size_t local_step = 0;
    bool running = true;

    for (Real t = 0;; ++step, ++local_step)
    {
        if (max_timesteps > 0 && local_step >= max_timesteps)
            break;

       // std::cout << "step " << step << std::endl;
        hpx::future<triple<Real> > local_max_velocity =
           part.do_timestep(dt);

        // if this is the root locality gather all remote residuals and sum up
        if (hpx::get_locality_id() == 0)
        {
            hpx::future<std::vector<triple<Real> > >
            max_velocities =
                hpx::lcos::gather_here(velocity_basename,
                                        std::move(local_max_velocity),
                                        num_localities, step);

            max_velocities.then(
                hpx::util::unwrapped(
                    [=, &t](std::vector<triple<Real> > local_max_uvws)
                    {
                        triple<Real> global_max_uvw(0);

                        for (auto& max_uvw : local_max_uvws)
                        {
                            global_max_uvw.x =
                                (max_uvw.x > global_max_uvw.x
                                    ? max_uvw.x : global_max_uvw.x);

                            global_max_uvw.y =
                                (max_uvw.y > global_max_uvw.y
                                    ? max_uvw.y : global_max_uvw.y);

                            global_max_uvw.z =
                                (max_uvw.z > global_max_uvw.z
                                    ? max_uvw.z : global_max_uvw.z);
                        }

                        Real new_dt =
                            std::min(re / 2. * 1. / (1. / std::pow(dx, 2)
                                        + 1. / std::pow(dy, 2)
                                        + 1. / std::pow(dz, 2))
                                    ,
                                    std::min(dx / global_max_uvw.x,
                                            std::min(dy / global_max_uvw.y, dz / global_max_uvw.z))
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
