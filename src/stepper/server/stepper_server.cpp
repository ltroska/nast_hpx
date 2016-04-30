#include <hpx/include/iostreams.hpp>
#include <hpx/lcos/gather.hpp>
#include <hpx/lcos/broadcast.hpp>

#include <cmath>

#include "stepper_server.hpp"

#include "io/vtk_writer.hpp"
#include "util/helpers.hpp"

typedef stepper::server::stepper_server stepper_component;
typedef hpx::components::component<stepper_component> stepper_server_type;

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_COMPONENT(stepper_server_type, stepper_component);
HPX_REGISTER_ACTION(stepper::server::stepper_server::setup_action,
    stepper_server_setup_action);

HPX_REGISTER_GATHER(RealType, stepper_server_space_gatherer);

namespace stepper
{
namespace server
{

stepper_server::stepper_server(uint nl)
: num_localities(nl)
{    
}

void stepper_server::setup(io::config&& cfg)
{
    // special case for two localities, we want a square configuration
    // of localities
    if (num_localities == 2)
    {
        num_localities_x = 2;
        num_localities_y = 1;
    }
    else
    {
        num_localities_x = static_cast<uint> (sqrt(num_localities));
        num_localities_y = num_localities_x;
    }

    c = std::move(cfg);

    // initialize things
    initialize_parameters();
    initialize_grids();
    initialize_communication();

    std::cout << cfg << std::endl;

    hpx::cout
        << "stepper on " << hpx::get_locality_id()
        << " with " << params.num_partitions_x - 2
        << "x" << params.num_partitions_y - 2 << " partitions, "
        << params.num_cells_per_partition_x
        << "x" << params.num_cells_per_partition_y 
        << " cells each, dx=" << params.dx 
        << " dy=" << params.dx
        << hpx::endl << hpx::flush;
    
    if (hpx::get_locality_id() == 0)
#ifdef CUSTOM_GRAIN_SIZE
        std::cout << "MODE: CUSTOM_GRAIN_SIZE" << std::endl;
#else
        std::cout << "MODE: WITH_FOR_EACH" << std::endl;
#endif

    communicate_uv_grid(0);

    // print out initial grids
    if ((c.output_skip_size != 0 || c.delta_vec != 0) && c.vtk)
        write_vtk(0);

    // start timestepping
    if (hpx::get_locality_id() == 0)
        do_work();
}

void stepper_server::initialize_parameters()
{
    params.num_cells_per_partition_x = c.i_res;
    params.num_cells_per_partition_y = c.j_res;

    // we want an even distribution of cells on each locality/partition
    params.num_partitions_x =
        ((c.i_max + 2) / num_localities_x) / c.i_res + 2;
    params.num_partitions_y =
        ((c.j_max + 2) / num_localities_y) / c.j_res + 2;

    if ((params.num_partitions_x - 2) * num_localities_x * c.i_res
            != c.i_max + 2
        ||
        (params.num_partitions_y - 2) * num_localities_y * c.j_res
            != c.j_max + 2)
    {
        std::cerr
            << "i_res and/or j_res doesn't fit the dimensions of the problem." 
            << std::endl;
        exit(0);
    }

    params.re = c.re;
    params.pr = c.pr;
    params.alpha = c.alpha;
    params.omega = c.omega;
    params.dx = c.x_length / c.i_max;
    params.dy = c.y_length / c.j_max;

    params.i_max = c.i_max;
    params.j_max = c.j_max;

    t = 0;
    next_write = 0;
    out_iter = 1;
}

void stepper_server::initialize_grids()
{
    index_grid.resize(params.num_partitions_x * params.num_partitions_y);

    uv_grid.resize(params.num_partitions_x * params.num_partitions_y);
    uv_temp_grid.resize(params.num_partitions_x * params.num_partitions_y);

    fg_grid.resize(params.num_partitions_x * params.num_partitions_y);
    fg_temp_grid.resize(params.num_partitions_x * params.num_partitions_y);

    p_grid.resize(params.num_partitions_x * params.num_partitions_y);
    p_temp_grid.resize(params.num_partitions_x * params.num_partitions_y);

    rhs_grid.resize(params.num_partitions_x * params.num_partitions_y);

    temperature_grid.resize(params.num_partitions_x * params.num_partitions_y);
    temperature_temp_grid.resize(
        params.num_partitions_x * params.num_partitions_y);

    stream_grid.resize(params.num_partitions_x * params.num_partitions_y);
    vorticity_grid.resize(params.num_partitions_x * params.num_partitions_y);
    heat_grid.resize(params.num_partitions_x * params.num_partitions_y);

    flag_grid.resize(params.num_partitions_x * params.num_partitions_y);

    scalar_dummy =
        scalar_partition(hpx::find_here(), params.num_cells_per_partition_x,
                            params.num_cells_per_partition_y);
    vector_dummy =
        vector_partition(hpx::find_here(), params.num_cells_per_partition_x,
                            params.num_cells_per_partition_y);

    for (uint l = 0; l < params.num_partitions_y; l++)
        for (uint k = 0; k < params.num_partitions_x; k++)
        {
            
            // fill grid with the global indices of bottom left cell for each
            // partition
            index_grid[get_index(k, l)] =
                std::pair<RealType, RealType>
                (
                    (hpx::get_locality_id() % num_localities_x)
                        * (params.num_partitions_x - 2)
                        * params.num_cells_per_partition_x
                        + (k - 1) * params.num_cells_per_partition_x
                ,
                    (hpx::get_locality_id() / num_localities_x)
                    * (params.num_partitions_y - 2)
                    * params.num_cells_per_partition_y
                    + (l - 1) * params.num_cells_per_partition_y
                );

            fg_grid[get_index(k, l)] =
                vector_partition(hpx::find_here(),
                    params.num_cells_per_partition_x,
                    params.num_cells_per_partition_y);
            
            p_grid[get_index(k, l)] =
                scalar_partition(hpx::find_here(),
                    params.num_cells_per_partition_x,
                    params.num_cells_per_partition_y,
                    0);
            
            rhs_grid[get_index(k, l)] =                
                scalar_partition(hpx::find_here(),
                params.num_cells_per_partition_x,
                params.num_cells_per_partition_y);
            
            temperature_grid[get_index(k, l)] =
                scalar_partition(hpx::find_here(),
                params.num_cells_per_partition_x,
                params.num_cells_per_partition_y, c.ti);
            
            stream_grid[get_index(k, l)] =
                scalar_partition(hpx::find_here(),
                params.num_cells_per_partition_x,
                params.num_cells_per_partition_y);
            
            vorticity_grid[get_index(k, l)] =
                scalar_partition(hpx::find_here(),
                params.num_cells_per_partition_x,
                params.num_cells_per_partition_y);
            
            heat_grid[get_index(k, l)] =
                scalar_partition(hpx::find_here(),
                params.num_cells_per_partition_x,
                params.num_cells_per_partition_y);

            std::vector<std::bitset<5> > local_flags(
                params.num_cells_per_partition_x
                * params.num_cells_per_partition_y);

            // get cell types from the types read from the csv file
            for (uint j = 0; j < params.num_cells_per_partition_y; j++)
                for (uint i = 0; i < params.num_cells_per_partition_x; i++)
                {
                    local_flags[j * params.num_cells_per_partition_x + i] =
                        c.flag_grid[
                            (params.j_max + 2 - 1 
                            - index_grid[get_index(k, l)].second - j)
                        *(params.i_max + 2)
                        + (index_grid[get_index(k, l)].first + i)
                        ];
                }

            flag_grid[get_index(k, l)] = std::move(local_flags);

            // populate velocities with initial data (if supplied)
            if (k > 0 && k < params.num_partitions_x - 1 && l > 0
                && l < params.num_partitions_y - 1)
            {
                if (c.with_initial_uv_grid)
                {
                    vector_data curr_data(params.num_cells_per_partition_x,
                        params.num_cells_per_partition_y);

                    for (uint j = 0; j < params.num_cells_per_partition_y; j++)
                        for (uint i = 0; i < params.num_cells_per_partition_x;
                                i++)
                        {
                            vector_cell& curr_cell =
                                curr_data.get_cell_ref(i, j);
                            
                            curr_cell.first =
                                c.initial_uv_grid[(params.j_max + 2 - 1
                                - index_grid[get_index(k, l)].second - j)
                                *(params.i_max + 2)
                                + (index_grid[get_index(k, l)].first + i)]
                                    .first;
                            
                            curr_cell.second =
                                c.initial_uv_grid[(params.j_max + 2 - 1
                                - index_grid[get_index(k, l)].second - j)
                                *(params.i_max + 2)
                                + (index_grid[get_index(k, l)].first + i)]
                                    .second;

                        }

                    uv_grid[get_index(k, l)] =
                        vector_partition(hpx::find_here(), curr_data);
                }
                else
                    uv_grid[get_index(k, l)] =
                        vector_partition(hpx::find_here(),
                            params.num_cells_per_partition_x,
                            params.num_cells_per_partition_y);
            }
            else
                uv_grid[get_index(k, l)] = vector_dummy;
        }
}

void stepper_server::initialize_communication()
{
    // find if locality has neighbor in each direction    
    has_neighbor[LEFT] = !(hpx::get_locality_id() % num_localities_x == 0);
    
    has_neighbor[RIGHT] =
        !(hpx::get_locality_id() % num_localities_x == num_localities_x - 1);
    
    has_neighbor[BOTTOM] = !(hpx::get_locality_id() / num_localities_x == 0);
    
    has_neighbor[TOP] =
        !(hpx::get_locality_id() / num_localities_x == num_localities_y - 1);
    
    has_neighbor[BOTTOM_LEFT] = (has_neighbor[BOTTOM] && has_neighbor[LEFT]);
    
    has_neighbor[BOTTOM_RIGHT] = (has_neighbor[BOTTOM] && has_neighbor[RIGHT]);
    
    has_neighbor[TOP_LEFT] = (has_neighbor[TOP] && has_neighbor[LEFT]);
    
    has_neighbor[TOP_RIGHT] = (has_neighbor[TOP] && has_neighbor[RIGHT]);

    
    for (int directionInt = LEFT; directionInt != NUM_DIRECTIONS;
                                                                directionInt++)
    {
        neighbor_steppers_[directionInt] =
            hpx::find_from_basename(stepper_basename,
                                    get_neighbor_id(hpx::get_locality_id(),
                                        static_cast<direction> (directionInt),
                                        num_localities));
    }

    // get the gids of remote steppers
    for (uint loc = 0; loc < num_localities; loc++)
        localities.push_back(
            hpx::find_from_basename(stepper_basename, loc).get());
}

void stepper_server::do_work()
{
    std::pair<RealType, RealType> max_uv(0, 0);
    RealType dt = c.dt;

    // start timestepping
    for (uint step = 1; t + dt < c.t_end; step++)
    {
        // tell all steppers to execute a timestep, returning the local maximal
        // velocities
        hpx::future<std::vector<std::pair<RealType, RealType> > > max_uvs =
            hpx::lcos::broadcast<do_timestep_action> (localities, step, dt);

        // when each locality finished, compute the global maximal velocity
        hpx::future<std::pair<RealType, RealType> > max_uv_fut =
            max_uvs.then(
            [](hpx::future<std::vector<std::pair<RealType, RealType> > > fut)
                -> std::pair<RealType, RealType>
            {
                auto local_max_uvs = fut.get();

                std::pair<RealType, RealType> res(0, 0);

                for (auto max_uv : local_max_uvs)
                {
                    res.first =
                        (max_uv.first > res.first ? max_uv.first : res.first);
                    
                    res.second =
                        (max_uv.second > res.second
                            ? max_uv.second : res.second);
                }

                return res;
            });

        max_uv = max_uv_fut.get();

        // CFL conditions
        RealType tmp =
            std::min(c.re / 2. * 1. / (1. / std::pow(params.dx, 2)
                        + 1. / std::pow(params.dy, 2))
                    ,
                    std::min(params.dx / max_uv.first,
                            params.dy / max_uv.second)
            );

        // special case for temperature driven flow
        if (c.pr)
            tmp = std::min(tmp, 
                            (c.re * c.pr) / 2. * 1.
                            / (1. / std::pow(params.dx, 2)
                            + 1. / std::pow(params.dy, 2)));

        dt = c.tau * tmp;
    }
}

std::pair<RealType, RealType> stepper_server::do_timestep(
    uint step, RealType dt)
{
    hpx::util::high_resolution_timer t1;
        
    // set the boundary and obstacle data for each partition
    for (uint l = 0; l < params.num_partitions_y; l++)
        for (uint k = 0; k < params.num_partitions_x; k++)
        {
            if (k != 0 && k != params.num_partitions_x - 1
                && l != 0 && l != params.num_partitions_y - 1)
            {
                uv_temp_grid[get_index(k, l)] =
                    hpx::dataflow(
                        hpx::launch::async,
                        &strategy::set_velocity_for_boundary_and_obstacles,
                        uv_grid[get_index(k, l)],
                        uv_grid[get_index(k - 1, l)],
                        uv_grid[get_index(k + 1, l)],
                        uv_grid[get_index(k, l - 1)],
                        uv_grid[get_index(k, l + 1)],
                        flag_grid[get_index(k, l)],
                        c.data_type,
                        c.u_bnd,
                        c.v_bnd
                    );

                if (c.pr != 0)
                    temperature_temp_grid[get_index(k, l)] =
                    hpx::dataflow(
                        hpx::launch::async,
                        &strategy::set_temperature_for_boundary_and_obstacles,
                        temperature_grid[get_index(k, l)],
                        temperature_grid[get_index(k - 1, l)],
                        temperature_grid[get_index(k + 1, l)],
                        temperature_grid[get_index(k, l - 1)],
                        temperature_grid[get_index(k, l + 1)],
                        flag_grid[get_index(k, l)],
                        c.temp_data_type,
                        c.temp_bnd,
                        index_grid[get_index(k, l)].first,
                        index_grid[get_index(k, l)].second,
                        params.dx,
                        params.dy
                    );

            }
            else
            {
                uv_temp_grid[get_index(k, l)] = uv_grid[get_index(k, l)];
                if (c.pr != 0)
                    temperature_temp_grid[get_index(k, l)] =
                        temperature_grid[get_index(k, l)];
            }
        }

    // note: this does not block (it is essentially moving on a future)
    uv_grid = uv_temp_grid;
  
    if (c.pr != 0)
        temperature_grid = temperature_temp_grid;

    communicate_uv_grid(step);

    // compute temperature on fluid cells if we have a temperature driven flow
    if (c.pr != 0)
    {
        for (uint l = 0; l < params.num_partitions_y; l++)
            for (uint k = 0; k < params.num_partitions_x; k++)
            {
                if (k != 0 && k != params.num_partitions_x - 1 && l != 0
                    && l != params.num_partitions_y - 1)
                {
                    temperature_temp_grid[get_index(k, l)] =
                        hpx::dataflow(
                            hpx::launch::async,
                            &strategy::compute_temperature_on_fluid_cells,
                            temperature_grid[get_index(k, l)],
                            temperature_grid[get_index(k - 1, l)],
                            temperature_grid[get_index(k + 1, l)],
                            temperature_grid[get_index(k, l - 1)],
                            temperature_grid[get_index(k, l + 1)],
                            uv_grid[get_index(k, l)],
                            uv_grid[get_index(k - 1, l)],
                            uv_grid[get_index(k, l - 1)],
                            flag_grid[get_index(k, l)],
                            params.re, c.pr, params.dx, params.dy, dt, c.alpha
                        );
                }
                else
                    temperature_temp_grid[get_index(k, l)] =
                        temperature_grid[get_index(k, l)];
            }


        temperature_grid = temperature_temp_grid;
    }

    // compute fg on fluid cells
    for (uint l = 1; l < params.num_partitions_y - 1; l++)
        for (uint k = 1; k < params.num_partitions_x - 1; k++)
        {
            fg_grid[get_index(k, l)]
                = hpx::dataflow(
                    hpx::launch::async,
                    &strategy::compute_fg_on_fluid_cells,
                    uv_grid[get_index(k, l)],
                    uv_grid[get_index(k - 1, l)],
                    uv_grid[get_index(k + 1, l)],
                    uv_grid[get_index(k, l - 1)],
                    uv_grid[get_index(k, l + 1)],
                    uv_grid[get_index(k + 1, l - 1)],
                    uv_grid[get_index(k - 1, l + 1)],
                    temperature_grid[get_index(k, l)],
                    temperature_grid[get_index(k + 1, l)],
                    temperature_grid[get_index(k, l + 1)],
                    flag_grid[get_index(k, l)],
                    c.re, c.gx, c.gy, c.beta, params.dx, params.dy, dt, c.alpha
                );
        }

   // communicate_fg_grid(step);

    // compute the right hand side for the Poisson equation of the pressure
    for (uint l = 1; l < params.num_partitions_y - 1; l++)
        for (uint k = 1; k < params.num_partitions_x - 1; k++)
            rhs_grid[get_index(k, l)] =
                hpx::dataflow(
                    hpx::launch::async,
                    &strategy::compute_right_hand_side_on_fluid_cells,
                    fg_grid[get_index(k, l)],
                    fg_grid[get_index(k - 1, l)],
                    fg_grid[get_index(k, l - 1)],
                    flag_grid[get_index(k, l)],
                    params.dx,
                    params.dy,
                    dt
                );
   
    /*print_grid(uv_grid, "uv");
    print_grid(temperature_grid, "temp");
    print_grid(fg_grid, "fg");
    print_grid(rhs_grid, "rhs");*/
    
    RealType t1_elapsed = t1.elapsed();

    hpx::util::high_resolution_timer t2;
    
    uint iter = 0;
    RealType res = 0;
    do
    {
        // if we were compiled with the hpx::parallel algorithm strategy,
        // we need to split setting the boundary values and doing the actual
        // SOR loop into two steps, otherwise we can do it in one
        
#ifdef CUSTOM_GRAIN_SIZE                 
        for (uint l = 0; l < params.num_partitions_y; l++)
            for (uint k = 0; k < params.num_partitions_x; k++)
            {
                if (k != 0 && k != params.num_partitions_x - 1 && l != 0
                    && l != params.num_partitions_y - 1)
                {
                    p_temp_grid[get_index(k, l)] =
                        hpx::dataflow(
                            hpx::launch::async,
                            &strategy::sor_cycle,
                            p_grid[get_index(k, l)],
                            p_temp_grid[get_index(k - 1, l)],
                            p_grid[get_index(k + 1, l)],
                            p_temp_grid[get_index(k, l - 1)],
                            p_grid[get_index(k, l + 1)],
                            rhs_grid[get_index(k, l)],
                            flag_grid[get_index(k, l)],
                            c.omega, params.dx, params.dy
                        );
                }
                else
                    p_temp_grid[get_index(k, l)] = p_grid[get_index(k, l)];
            }
        
        p_grid = p_temp_grid;        
#else
        for (uint l = 0; l < params.num_partitions_y; l++)
            for (uint k = 0; k < params.num_partitions_x; k++)
            {
                if (k != 0 && k != params.num_partitions_x - 1 && l != 0
                    && l != params.num_partitions_y - 1)
                {
                    p_temp_grid[get_index(k, l)] =
                        hpx::dataflow(
                            hpx::launch::async,
                            &strategy::set_pressure_for_boundary_and_obstacles,
                            p_grid[get_index(k, l)],
                            p_grid[get_index(k - 1, l)],
                            p_grid[get_index(k + 1, l)],
                            p_grid[get_index(k, l - 1)],
                            p_grid[get_index(k, l + 1)],
                            flag_grid[get_index(k, l)]
                        );
                }
                else
                    p_temp_grid[get_index(k, l)] = p_grid[get_index(k, l)];
            }
        
        for (uint l = 0; l < params.num_partitions_y; l++)
            for (uint k = 0; k < params.num_partitions_x; k++)
            {
                if (k != 0 && k != params.num_partitions_x - 1 && l != 0
                    && l != params.num_partitions_y - 1)
                {
                    p_grid[get_index(k, l)] =
                        hpx::dataflow(
                            hpx::launch::async,
                            &strategy::sor_cycle,
                            p_temp_grid[get_index(k, l)],
                            p_grid[get_index(k - 1, l)],
                            p_temp_grid[get_index(k + 1, l)],
                            p_grid[get_index(k, l - 1)],
                            p_temp_grid[get_index(k, l + 1)],
                            rhs_grid[get_index(k, l)],
                            flag_grid[get_index(k, l)],
                            c.omega, params.dx, params.dy
                        );
                }
                else
                    p_grid[get_index(k, l)] = p_temp_grid[get_index(k, l)];
            }
#endif
        
        communicate_p_grid(step * c.iter_max + iter);
          
        std::vector<hpx::future<RealType> > residuals;

        for (uint l = 1; l < params.num_partitions_y - 1; l++)
            for (uint k = 1; k < params.num_partitions_x - 1; k++)
                residuals.push_back(
                        hpx::dataflow(
                            hpx::launch::async,
                            &strategy::compute_residual,
                            p_grid[get_index(k, l)],
                            p_grid[get_index(k - 1, l)],
                            p_grid[get_index(k + 1, l)],
                            p_grid[get_index(k, l - 1)],
                            p_grid[get_index(k, l + 1)],
                            rhs_grid[get_index(k, l)],
                            flag_grid[get_index(k, l)],
                            params.dx, params.dy
                        )
                );

        // sum up all local residuals
        hpx::future<RealType> residual_fut =
            hpx::when_all(residuals).then(
                hpx::util::unwrapped(
                    [](std::vector< hpx::future<RealType> > residuals)
                    -> RealType
                    {
                        RealType sum = 0;

                        for (auto& r : residuals)
                            sum += r.get();

                        return sum;
                    }
                )
            );

        iter++;
        
        // if this is the root locality gather all remote residuals and sum up
        if (hpx::get_locality_id() == 0)
        {
            hpx::future<std::vector<RealType> > local_residuals =
                hpx::lcos::gather_here(gather_basename,
                                        std::move(residual_fut),
                                        num_localities, step*c.iter_max + iter);

            hpx::future<RealType> residual =
                local_residuals.then(
                    [](hpx::future<std::vector<RealType>> local_residuals)
                        -> RealType
                    {
                        RealType result = 0;
                        std::vector<RealType> local_res = local_residuals.get();

                        for (RealType res : local_res)
                            result += res;

                        return result;
                    });

            res = residual.get() / c.num_fluid_cells;
            
            // decide if SOR should keep running or not
            hpx::lcos::broadcast_apply<set_keep_running_action>(
                localities, step * c.iter_max + iter,
                (iter < c.iter_max && res > c.eps_sq));
        }
        // if not root locality, send residual to root locality
        else
            hpx::lcos::gather_there(gather_basename, std::move(residual_fut),
                                        step*c.iter_max + iter).wait();
    }
    while (keep_running.receive(step * c.iter_max + iter).get());
    
   // print_grid(p_grid, "p");
    
    RealType t2_elapsed = t2.elapsed();

    hpx::util::high_resolution_timer t3;

    // update velocities and find local maximal velocity
    
    std::vector<hpx::future<std::pair<RealType, RealType> > > max_uvs;

    for (uint l = 1; l < params.num_partitions_y - 1; l++)
        for (uint k = 1; k < params.num_partitions_x - 1; k++)
        {
            hpx::future<std::pair<vector_partition,
                                    std::pair<RealType, RealType> > >
                f = hpx::dataflow(
                        hpx::launch::async,
                        &strategy::update_velocities,
                        uv_grid[get_index(k, l)],
                        p_grid[get_index(k, l)],
                        p_grid[get_index(k + 1, l)],
                        p_grid[get_index(k, l + 1)],
                        fg_grid[get_index(k, l)],
                        flag_grid[get_index(k, l)],
                        params.dx, params.dy, dt
                    );

            // we need to split the returned future<pair> into two futures
            // to be able to use them independently
            hpx::lcos::local::promise<std::pair<RealType, RealType> > outer_p;

            max_uvs.emplace_back(outer_p.get_future());

            uv_grid[get_index(k, l)] = f.then(
                [outer_p = std::move(outer_p)](auto fut) mutable
                -> vector_partition
                {
                    std::pair<vector_partition,
                            std::pair<RealType, RealType> >
                        p = fut.get();

                    outer_p.set_value(p.second);

                    return p.first;
                });

        }

   // print_grid(uv_grid, "uv2");
    
    t += dt;

    // compute local maximal velocities
    hpx::future<std::pair<RealType, RealType> > max_uv =
        hpx::when_all(max_uvs).then(
        hpx::util::unwrapped(
        [](std::vector< hpx::future<std::pair<RealType, RealType> > > max_uvs)
        -> std::pair<RealType, RealType>
        {
            std::pair<RealType, RealType> max_uv(0, 0);

            for (auto& m_uv_f : max_uvs)
            {
                auto m_uv = m_uv_f.get();

                    max_uv.first =
                        (m_uv.first > max_uv.first ?
                        m_uv.first : max_uv.first);

                    max_uv.second =
                        (m_uv.second > max_uv.second ?
                        m_uv.second : max_uv.second);
            }

            return max_uv;
        })
    );
        
    auto max_ = max_uv.get();

    // print out local grid
    if ((c.output_skip_size != 0 && (step % c.output_skip_size == 0))
        || (c.delta_vec != 0 && next_write <= t))
    {
        if (c.delta_vec != 0)
            next_write = t + c.delta_vec;

        if (c.vtk)
            write_vtk(out_iter);

        if (hpx::get_locality_id() == 0)
            std::cout 
                << "t " << t << " | dt " << dt 
                << " | iterations: " << iter 
                << " | residual squared " << res
                << " | before SOR = " << t1_elapsed 
                << " | SOR = " << t2_elapsed 
                << " average = " << t2_elapsed / iter 
                << " | after SOR = " << t3.elapsed() 
                << " | max uv " << max_.first << " " << max_.second
                << std::endl;

        out_iter++;
    }

    return max_;
}

// ---------------------------- COMMUNICATION ------------------------------ //

void stepper_server::set_keep_running(uint iter, bool kr)
{
    keep_running.store_received(iter, std::move(kr));
}

void stepper_server::communicate_p_grid(uint iter)
{
    // send data from all four sides in corresponding direction
    for (uint l = 1; l < params.num_partitions_y - 1; l++)
    {
        send_p_to_neighbor(iter * params.num_partitions_y + l,
            p_grid[get_index(1, l)], LEFT);
        
        send_p_to_neighbor(iter * params.num_partitions_y + l,
            p_grid[get_index(params.num_partitions_x - 2, l)], RIGHT);
    }

    for (uint k = 1; k < params.num_partitions_x - 1; k++)
    {
        send_p_to_neighbor(iter * params.num_partitions_x + k,
            p_grid[get_index(k, 1)], BOTTOM);
        
        send_p_to_neighbor(iter * params.num_partitions_x + k,
            p_grid[get_index(k, params.num_partitions_y - 2)], TOP);
    }

    // receive data from localities neighboring all four sides
    for (uint l = 1; l < params.num_partitions_y - 1; l++)
    {
        p_grid[get_index(0, l)] =
            receive_p_from_neighbor(iter * params.num_partitions_y + l, LEFT);
        
        p_grid[get_index(params.num_partitions_x - 1, l)] =
            receive_p_from_neighbor(iter * params.num_partitions_y + l, RIGHT);
    }

    for (uint k = 1; k < params.num_partitions_x - 1; k++)
    {
        p_grid[get_index(k, 0)] =
            receive_p_from_neighbor(iter * params.num_partitions_x + k, BOTTOM);
        
        p_grid[get_index(k, params.num_partitions_y - 1)] =
            receive_p_from_neighbor(iter * params.num_partitions_x + k, TOP);
    }
}

void stepper_server::receive_p_action_(
    uint t, scalar_partition p, direction to_dir)
{
    //we need to negate the direction to put the data into the correct buffer
    p_recv_buffs_[NUM_DIRECTIONS - to_dir - 1].store_received(t, std::move(p));
}

void stepper_server::send_p_to_neighbor(
    uint t, scalar_partition p, direction dir)
{
    if (has_neighbor[dir])
    {
        receive_p_action act;
        hpx::async(act, neighbor_steppers_[dir].get(), t, p, dir);
    }
}

scalar_partition stepper_server::receive_p_from_neighbor(uint t, direction dir)
{
    if (has_neighbor[dir])
        return p_recv_buffs_[dir].receive(t);
    else
        // return dummy to fake neighboring locality
        return scalar_dummy;
}

void stepper_server::communicate_fg_grid(uint step)
{
    // send data from all four sides in corresponding direction
    for (uint l = 1; l < params.num_partitions_y - 1; l++)
        send_fg_to_neighbor(step * params.num_partitions_y + l,
            fg_grid[get_index(params.num_partitions_x - 2, l)], RIGHT);

    for (uint k = 1; k < params.num_partitions_x - 1; k++)
        send_fg_to_neighbor(step * params.num_partitions_x + k,
            fg_grid[get_index(k, params.num_partitions_y - 2)], TOP);

    // receive data from localities neighboring all four sides
    for (uint l = 1; l < params.num_partitions_y - 1; l++)
        fg_grid[get_index(0, l)] =
            receive_fg_from_neighbor(step * params.num_partitions_y + l, LEFT);

    for (uint k = 1; k < params.num_partitions_x - 1; k++)
        fg_grid[get_index(k, 0)] =
            receive_fg_from_neighbor(step * params.num_partitions_x + k, BOTTOM);
}

void stepper_server::receive_fg_action_(
    uint t, vector_partition fg, direction to_dir)
{
    //we need to negate the direction to put the data into the correct buffer
    fg_recv_buffs_[NUM_DIRECTIONS - to_dir - 1].store_received(t,
                                                                std::move(fg));
}

void stepper_server::send_fg_to_neighbor(
    uint t, vector_partition fg, direction dir)
{
    if (has_neighbor[dir])
    {
        receive_fg_action act;
        hpx::async(act, neighbor_steppers_[dir].get(), t, fg, dir);
    }
}

vector_partition stepper_server::receive_fg_from_neighbor(uint t, direction dir)
{
    if (has_neighbor[dir])
        return fg_recv_buffs_[dir].receive(t);
    else
        // return dummy to fake neighboring locality
        return vector_dummy;
}

uint stepper_server::get_index(uint k, uint l) const
{
    return l * params.num_partitions_x + k;
}

void stepper_server::communicate_uv_grid(uint step)
{
    // SEND
    for (uint l = 1; l < params.num_partitions_y - 1; l++)
    {
        if (l == 1)
            send_uv_to_neighbor(step,
                uv_grid[get_index(params.num_partitions_x - 2, l)],
                BOTTOM_RIGHT);

        if (l == params.num_partitions_y - 2)
            send_uv_to_neighbor(step, uv_grid[get_index(1, l)], TOP_LEFT);

        send_uv_to_neighbor(step * params.num_partitions_y + l,
            uv_grid[get_index(1, l)], LEFT);
        
        send_uv_to_neighbor(step * params.num_partitions_y + l,
            uv_grid[get_index(params.num_partitions_x - 2, l)], RIGHT);
    }

    for (uint k = 1; k < params.num_partitions_x - 1; k++)
    {
        if (k == 1)
            send_uv_to_neighbor(step, uv_grid[get_index(k, 1)], BOTTOM_LEFT);

        if (k == params.num_partitions_x - 2)
            send_uv_to_neighbor(step, uv_grid[get_index(k,
                params.num_partitions_y - 2)], TOP_RIGHT);

        send_uv_to_neighbor(step * params.num_partitions_x + k,
            uv_grid[get_index(k, 1)], BOTTOM);
        
        send_uv_to_neighbor(step * params.num_partitions_x + k,
            uv_grid[get_index(k, params.num_partitions_y - 2)], TOP);
    }

    // RECEIVE
    for (uint l = 1; l < params.num_partitions_y - 1; l++)
    {
        uv_grid[get_index(0, l)] =
            receive_uv_from_neighbor(step * params.num_partitions_y + l, LEFT);
        
        uv_grid[get_index(params.num_partitions_x - 1, l)] =
            receive_uv_from_neighbor(step * params.num_partitions_y + l, RIGHT);
    }

    for (uint k = 1; k < params.num_partitions_x - 1; k++)
    {
        uv_grid[get_index(k, 0)] =
            receive_uv_from_neighbor(step * params.num_partitions_x + k,
                                                                        BOTTOM);
        
        uv_grid[get_index(k, params.num_partitions_y - 1)] =
            receive_uv_from_neighbor(step * params.num_partitions_x + k, TOP);
    }

    uv_grid[get_index(0, params.num_partitions_y - 1)] =
        receive_uv_from_neighbor(step, TOP_LEFT);
    
    uv_grid[get_index(params.num_partitions_x - 1,
                        params.num_partitions_y - 1)]
        = receive_uv_from_neighbor(step, TOP_RIGHT);
    
    uv_grid[get_index(0, 0)] = receive_uv_from_neighbor(step, BOTTOM_LEFT);
    
    uv_grid[get_index(params.num_partitions_x - 1, 0)] =
        receive_uv_from_neighbor(step, BOTTOM_RIGHT);
}

void stepper_server::receive_uv_action_(
    uint t, vector_partition uv, direction to_dir)

{   
    //we need to negate the direction to put the data into the correct buffer
    uv_recv_buffs_[NUM_DIRECTIONS - to_dir - 1].store_received(t,
                                                                std::move(uv));
}

void stepper_server::send_uv_to_neighbor(
    uint t, vector_partition uv, direction dir)
{
    if (has_neighbor[dir])
    {
        receive_uv_action act;
        hpx::async(act, neighbor_steppers_[dir].get(), t, uv, dir);
    }
}

vector_partition stepper_server::receive_uv_from_neighbor(uint t, direction dir)
{
    if (has_neighbor[dir])
        return uv_recv_buffs_[dir].receive(t);
    else
        // return dummy to fake neighboring locality
        return vector_dummy;
}

template<typename T>
void stepper_server::print_grid(
    std::vector<grid::partition<T> > const& grid, const std::string message)
const
{
    std::shared_ptr<hpx::lcos::local::promise<int> > p =
        std::make_shared<hpx::lcos::local::promise<int> >();
    
    io::do_async_print(grid, message, params.num_partitions_x,
        params.num_partitions_y, params.num_cells_per_partition_x,
        params.num_cells_per_partition_y, p);
}

void stepper_server::write_vtk(uint step)
{
    // compute the visualization data
    for (uint l = 1; l < params.num_partitions_y - 1; l++)
        for (uint k = 1; k < params.num_partitions_x - 1; k++)
        {
            uint global_i = index_grid[get_index(k, l)].first;
            uint global_j = index_grid[get_index(k, l)].second;

            hpx::future<std::tuple<scalar_partition,
                                    scalar_partition,
                                    scalar_partition> >
                f = hpx::dataflow(
                        hpx::launch::async,
                        &strategy::compute_stream_vorticity_heat,
                        stream_grid[get_index(k, l)],
                        stream_grid[get_index(k, l - 1)],
                        vorticity_grid[get_index(k, l)],
                        heat_grid[get_index(k, l)],
                        heat_grid[get_index(k, l - 1)],
                        uv_grid[get_index(k, l)],
                        uv_grid[get_index(k + 1, l)],
                        uv_grid[get_index(k, l +1)],
                        temperature_grid[get_index(k, l)],
                        temperature_grid[get_index(k + 1, l)],
                        flag_grid[get_index(k, l)],
                        global_i, global_j, params.i_max,
                        params.j_max, c.re, c.pr, params.dx, params.dy
                    );
            
            hpx::lcos::local::promise<scalar_partition> outer_p1;
            hpx::lcos::local::promise<scalar_partition> outer_p2;
            
            stream_grid[get_index(k, l)] = f.then(
                [outer_p1 = std::move(outer_p1), outer_p2 = std::move(outer_p2)]
                (auto fut) mutable -> scalar_partition
                {
                    std::tuple<scalar_partition, scalar_partition,
                                scalar_partition>
                        p = fut.get();
                    
                    outer_p1.set_value(std::get<1>(p));
                    outer_p2.set_value(std::get<2>(p));
                    
                    return std::get<0>(p);
                }
            );
        }
        
    // write out data once computation finishes
    hpx::dataflow(
        hpx::launch::async,
        &io::write_vtk,
        p_grid, uv_grid, stream_grid, vorticity_grid, heat_grid,
        temperature_grid, flag_grid, params.dx, params.dy, step, params.i_max,
        params.j_max, params.num_partitions_x,
        params.num_partitions_y, params.num_cells_per_partition_x,
        params.num_cells_per_partition_y
    );
}

}//namespace server
}//namespace stepper
