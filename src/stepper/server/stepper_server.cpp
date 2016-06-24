#include <chrono>

//#include <hpx/include/iostreams.hpp>
#include <hpx/lcos/gather.hpp>
#include <hpx/lcos/broadcast.hpp>
#include <hpx/runtime/components/migrate_component.hpp>
#include <hpx/lcos/barrier.hpp>
#include <hpx/lcos/wait_all.hpp>
#include <hpx/parallel/algorithms/for_each.hpp>
#include <hpx/parallel/algorithms/transform_reduce.hpp>

#include <boost/range/irange.hpp>

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
    auto here = hpx::get_locality_id();   

    auto num_localities = hpx::get_num_localities_sync();
    
    hpx::lcos::barrier b;    
    if (here == 0)
    {
         
        b = std::move(hpx::lcos::barrier::create(hpx::find_here(), num_localities));
        hpx::agas::register_name_sync(barrier_basename, b.get_id());
    }
    else
    {
        hpx::id_type idb = hpx::agas::on_symbol_namespace_event(
                barrier_basename, hpx::agas::symbol_ns_bind, true).get();
        b = std::move(hpx::lcos::barrier(idb));
    }
        
    Real dt = cfg.initial_dt;
    Real dx = cfg.dx;
    Real dy = cfg.dy;
    Real re = cfg.re;
    Real pr = cfg.pr;
    Real tau = cfg.tau;
    Real t_end = cfg.t_end;

    std::vector<grid::partition> parts;
    for (std::size_t idy = 0; idy < cfg.num_local_partitions_y; ++idy)
        for (std::size_t idx = 0; idx < cfg.num_local_partitions_x; ++idx)
        {
            auto act_idx = (here % cfg.num_localities_x) * cfg.num_local_partitions_x + idx;
            auto act_idy = (here / cfg.num_localities_x) * cfg.num_local_partitions_y + idy;
            
            auto local_idx = idx;
            auto local_idy = idy;            
            
            auto rank = act_idy * cfg.num_local_partitions_x * cfg.num_localities_x + act_idx;
                        
            parts.emplace_back(
                grid::partition(hpx::find_here(), cfg, act_idx, act_idy, local_idx, local_idy, rank));
                
           // parts[idy * cfg.num_local_partitions_x + idx].init_sync();
        }
        
    auto rge = boost::irange(static_cast<std::size_t>(0), cfg.num_local_partitions);
        
    hpx::parallel::for_each(hpx::parallel::par, std::begin(rge), std::end(rge),
        [&](std::size_t p)
        {
            parts[p].init_sync();
        }
    );
    
    std::size_t step = 0;

    std::vector<hpx::naming::id_type> localities;
    for (uint loc = 0; loc < num_localities; loc++)
        localities.push_back(hpx::find_from_basename(stepper_basename, loc).get());

    if (hpx::get_locality_id() == 1)
    {
       // std::cout << "migrating " << std::endl;
       // if (pr == 0)
    //    hpx::components::migrate(parts[0], hpx::find_remote_localities()[0]).wait();
    //    parts[0].init_sync();
    }

    b.wait();
        
    bool running = true;
    for (Real t = 0;; step++)
    {
       // std::cout << "step " << step << std::endl;
       
        hpx::future<std::pair<Real, Real> > local_max_velocity =
            hpx::parallel::transform_reduce(hpx::parallel::par(hpx::parallel::task), std::begin(rge), std::end(rge),
                [&](std::size_t p) -> std::pair<Real, Real>
                {
                    return parts[p].do_timestep(dt).get();
                },
                std::make_pair(0., 0.),
                [](std::pair<Real, Real> a, std::pair<Real, Real> b) -> std::pair<Real, Real>
                {
                    return std::make_pair(a.first > b.first ? a.first : b.first,
                        a.second > b.second ? a.second : b.second);
                }
            );

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

        if (t >= t_end)
            break;
        t += dt;
        dt = dt_buffer.receive(step).get();

        if (hpx::get_locality_id() == 0)
            for (std::size_t loc = 0; loc < num_localities; ++ loc)
            {
                hpx::performance_counters::performance_counter count(
                    "/threads{locality#" + std::to_string(loc) + "/total}/idle-rate");
                std::cout << loc << " " << count.get_value<double>().get() << std::endl;
            }
    }
}

void stepper_server::do_timestep(Real dt)
{

}

void stepper_server::set_dt(uint step, Real dt)
{
    dt_buffer.store_received(step, std::move(dt));
}


}//namespace server
}//namespace stepper
}
