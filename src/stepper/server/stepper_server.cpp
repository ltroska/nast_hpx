#include "stepper_server.hpp"

#include "util/triple.hpp"

#include <chrono>

typedef nast_hpx::stepper::server::stepper_server stepper_component;
typedef hpx::components::component<stepper_component> stepper_server_type;

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_COMPONENT(stepper_server_type, stepper_component);
HPX_REGISTER_ACTION(nast_hpx::stepper::server::stepper_server::setup_action,
    stepper_server_setup_action);
HPX_REGISTER_ACTION(nast_hpx::stepper::server::stepper_server::run_action,
    stepper_server_run_action);

typedef nast_hpx::triple<double> vec3;
HPX_REGISTER_GATHER(vec3, stepper_server_velocity_gather);

namespace nast_hpx { namespace stepper { namespace server {

stepper_server::stepper_server(uint nl)
: num_localities(nl)
{}

void stepper_server::setup(io::config const& cfg)
{
    rank = hpx::get_locality_id();
    num_localities = hpx::get_initial_num_localities();

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

    std::vector<hpx::future<hpx::id_type > > steps =
        hpx::find_all_from_basename(stepper_basename, num_localities);

    localities = hpx::when_all(steps).then(hpx::util::unwrapped2(
                [](std::vector<hpx::id_type>&& ids) -> std::vector<hpx::id_type>
                { return ids;})
            ).get();
}

void stepper_server::run()
{
    part.init_sync();

    double dt = init_dt;

    std::size_t local_step = 0;

    for (double t = 0;; ++step, ++local_step)
    {
        if (max_timesteps > 0 && local_step >= max_timesteps)
            break;

        hpx::future<triple<double> > local_max_velocity =
           part.do_timestep(dt);

        // if this is the root locality gather all remote residuals and sum up
        if (hpx::get_locality_id() == 0)
        {
            hpx::future<std::vector<triple<double> > >
            max_velocities =
                hpx::lcos::gather_here(velocity_basename,
                                        std::move(local_max_velocity),
                                        num_localities, step);

            max_velocities.then(
                hpx::util::unwrapped(
                    [=, &t](std::vector<triple<double> > local_max_uvws)
                    {
                        triple<double> global_max_uvw(0);

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

                        double new_dt =
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

        if (t >= t_end)
            break;
        t += dt;
        dt = dt_buffer.receive(step).get();
    }
}

void stepper_server::set_dt(uint step, double dt)
{
    dt_buffer.store_received(step, std::move(dt));
}


}//namespace server
}//namespace stepper
}
