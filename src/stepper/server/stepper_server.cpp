#include <hpx/include/iostreams.hpp>
#include <hpx/lcos/gather.hpp>

#include <cmath>

#include "stepper_server.hpp"
#include "computation/custom_grain_size.hpp"

#include "io/vtk_writer.hpp"

typedef stepper::server::stepper_server stepper_component;
typedef hpx::components::component<stepper_component> stepper_server_type;

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_COMPONENT(stepper_server_type, stepper_component);
HPX_REGISTER_ACTION(stepper::server::stepper_server::setup_action, stepper_server_setup_action);

HPX_REGISTER_GATHER(RealType, stepper_server_space_gatherer)

namespace stepper { namespace server {

stepper_server::stepper_server(uint nl)
    : num_localities(nl)
{}

void stepper_server::setup(io::config cfg)
{
    c = cfg;
    initialize_parameters();
    initialize_grids();
    initialize_communication();

    strategy = new computation::custom_grain_size(index_grid, params);

    std::cout << "stepper on " << hpx::get_locality_id() << " with " <<  params.num_partitions_x-2 << "x"<< params.num_partitions_y-2 << " partitions, "
             << params.num_cells_per_partition_x << "x" <<  params.num_cells_per_partition_y << " cells each, " << "dx=" << params.dx << " dy=" << params.dx
             << std::endl;

    do_work();

}

void stepper_server::initialize_parameters()
{
    if (num_localities == 2)
    {
        num_localities_x = 2;
        num_localities_y = 1;
    }
    else
    {
        num_localities_x = static_cast<uint>(sqrt(num_localities));
        num_localities_y = num_localities_x;
    }

    params.num_cells_per_partition_x = c.i_res;
    params.num_cells_per_partition_y = c.j_res;

    params.num_partitions_x = ((c.i_max + 2) / num_localities_x) / c.i_res + 2;
    params.num_partitions_y = ((c.j_max + 2) / num_localities_y) / c.j_res + 2;

    params.re = c.re;
    params.alpha = c.alpha;
    params.omega = c.omega;
    params.dx = c.x_length / c.i_max;
    params.dy = c.y_length / c.j_max;

    params.i_max = c.i_max;
    params.j_max = c.j_max;
}

void stepper_server::initialize_grids()
{
    index_grid.resize(params.num_partitions_x * params.num_partitions_y);
    uv_grid.resize(params.num_partitions_x * params.num_partitions_y);
    fg_grid.resize(params.num_partitions_x * params.num_partitions_y);
    p_grid.resize(params.num_partitions_x * params.num_partitions_y);
    rhs_grid.resize(params.num_partitions_x * params.num_partitions_y);

    for (uint l = 0; l < params.num_partitions_y; l++)
        for (uint k = 0; k < params.num_partitions_x; k++)
        {
            index_grid[get_index(k, l)] =
                std::pair<RealType, RealType>
                    (
                    (hpx::get_locality_id() % num_localities_x) * (params.num_partitions_x - 2) * params.num_cells_per_partition_x + (k - 1) * params.num_cells_per_partition_x,
                    (hpx::get_locality_id() / num_localities_x) * (params.num_partitions_y - 2) * params.num_cells_per_partition_y + (l - 1) * params.num_cells_per_partition_y
                    );

            uv_grid[get_index(k, l)] = vector_partition(hpx::find_here(), params.num_cells_per_partition_x, params.num_cells_per_partition_y);
            fg_grid[get_index(k, l)] = vector_partition(hpx::find_here(), params.num_cells_per_partition_x, params.num_cells_per_partition_y);
            p_grid[get_index(k, l)] = scalar_partition(hpx::find_here(), params.num_cells_per_partition_x, params.num_cells_per_partition_y);
            rhs_grid[get_index(k, l)] = scalar_partition(hpx::find_here(), params.num_cells_per_partition_x, params.num_cells_per_partition_y);
        }

    scalar_dummy = scalar_partition(hpx::find_here(), 1, 1);
    vector_dummy = vector_partition(hpx::find_here(), 1, 1);
}

void stepper_server::initialize_communication()
{
    has_neighbor[LEFT] = !(hpx::get_locality_id() % num_localities_x == 0);
    has_neighbor[RIGHT] = !(hpx::get_locality_id() % num_localities_x == num_localities_x - 1);
    has_neighbor[BOTTOM] = !(hpx::get_locality_id() / num_localities_x == 0);
    has_neighbor[TOP] = !(hpx::get_locality_id() / num_localities_x == num_localities_y - 1);
    has_neighbor[BOTTOM_LEFT] = (has_neighbor[BOTTOM] && has_neighbor[LEFT]);
    has_neighbor[BOTTOM_RIGHT] = (has_neighbor[BOTTOM] && has_neighbor[RIGHT]);
    has_neighbor[TOP_LEFT] = (has_neighbor[TOP] && has_neighbor[LEFT]);
    has_neighbor[TOP_RIGHT] = (has_neighbor[TOP] && has_neighbor[RIGHT]);

    for (int directionInt = LEFT; directionInt != NUM_DIRECTIONS; directionInt++)
    {
        neighbor_steppers_[directionInt] = hpx::find_from_basename(stepper_basename,
                                                                   get_neighbor_id(hpx::get_locality_id(), static_cast<direction>(directionInt), num_localities));
        std::cout << get_neighbor_id(hpx::get_locality_id(), static_cast<direction>(directionInt), num_localities) << " ";
    }

    std::cout << std::endl;

}

void stepper_server::do_work()
{
    std::pair<RealType, RealType> max_uv(2, 0);

    RealType t = 0, dt = 0;
    for (uint step = 0; t < c.t_end; step++)
    {
        dt = compute_new_dt(max_uv);
        t += dt;

     //   std::cout << "step: " << step << " t: " << t << " dt: " << dt;
        max_uv = do_timestep(step, dt).get();
     //   std:: cout << " max_uv: " << max_uv.first << " " << max_uv.second << std::endl;

       // print_grid(p_grid);
      // if (step % c.output_skip_size == 0)
        //    write_vtk(step);
        break;
    }
}

hpx::future<std::pair<RealType, RealType> > stepper_server::do_timestep(uint step, RealType dt)
{
    strategy->set_velocity_on_boundary(uv_grid);

    communicate_uv_grid(step);

    strategy->compute_fg(fg_grid, uv_grid, dt);

    communicate_fg_grid(step);

    strategy->compute_rhs(rhs_grid, fg_grid, dt);

    RealType residual;
    uint iter;
    for (iter = 0; iter < c.iter_max; iter++)
    {
        std::cout << "iteration: " << iter << std::endl;
        strategy->set_pressure_on_boundary(p_grid);

        hpx::future<RealType> residual_fut = strategy->sor_cycle(p_grid, rhs_grid);
        communicate_p_grid(iter);
      //  std::cout << "iteration: " << iter << " | local residual: " << residual << std::endl<<std::flush;
        if (hpx::get_locality_id() == 0)
        {
        hpx::future<std::vector<RealType> > overall_result =
            hpx::lcos::gather_here(gather_basename, std::move(residual_fut), num_localities);

            std::vector<RealType> sol = overall_result.get();

            std::cout << sol.size();


            for (int i = 0; i < sol.size(); i++)
            {
                std::cout << "got residual " << sol[i]  << " from " << i << std::endl;
            }
        }
        else
        {
            hpx::lcos::gather_there(gather_basename, std::move(residual_fut)).wait();
        }

       // if (residual < c.eps_sq)
         //   break;
    }

   // std::cout << " iterations: " << iter << " residual: " << residual;

    return strategy->update_velocities(uv_grid, fg_grid, p_grid, dt);
}

// ---------------------------------------- COMMUNICATION ---------------------------------------- //
void stepper_server::communicate_p_grid(uint iter)
{
    //SEND
    for (uint l = 1; l < params.num_partitions_y - 1; l++)
    {
        send_p_to_neighbor(iter * params.num_partitions_y + l, p_grid[get_index(1, l)], LEFT);
        send_p_to_neighbor(iter * params.num_partitions_y + l, p_grid[get_index(params.num_partitions_x - 2, l)], RIGHT);
    }

    for (uint k = 1; k < params.num_partitions_x - 1; k++)
    {
        send_p_to_neighbor(iter * params.num_partitions_x + k, p_grid[get_index(k, 1)], BOTTOM);
        send_p_to_neighbor(iter * params.num_partitions_x + k, p_grid[get_index(k, params.num_partitions_y - 2)], TOP);
    }

    //RECEIVE
    for(uint l = 1; l < params.num_partitions_y - 1; l++)
    {
        p_grid[get_index(0, l)] = receive_p_from_neighbor(iter * params.num_partitions_y+l, LEFT);
        p_grid[get_index(params.num_partitions_x - 1, l)] = receive_p_from_neighbor(iter * params.num_partitions_y + l, RIGHT);
    }

    for(uint k = 1; k < params.num_partitions_x - 1; k++)
    {
        p_grid[get_index(k, 0)] = receive_p_from_neighbor(iter * params.num_partitions_x + k, BOTTOM);
        p_grid[get_index(k, params.num_partitions_y - 1)] = receive_p_from_neighbor(iter * params.num_partitions_x + k, TOP);
    }
}

void stepper_server::receive_p_action_(uint t, scalar_partition p, direction to_dir)
{
    //direction is now opposite.
    p_recv_buffs_[NUM_DIRECTIONS - to_dir - 1].store_received(t, std::move(p));
}

void stepper_server::send_p_to_neighbor(uint t, scalar_partition p, direction dir)
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
        return scalar_dummy;
}

void stepper_server::communicate_fg_grid(uint step)
{
    //SEND
     for (uint l = 1; l < params.num_partitions_y - 1; l++)
        send_fg_to_neighbor(step * params.num_partitions_y + l, fg_grid[get_index(params.num_partitions_x - 2, l)], RIGHT);

    for (uint k = 1; k < params.num_partitions_x - 1; k++)
        send_fg_to_neighbor(step * params.num_partitions_x + k, fg_grid[get_index(k, params.num_partitions_y - 2)], TOP);

    //RECEIVE
    for(uint l = 1; l < params.num_partitions_y - 1; l++)
        fg_grid[get_index(0, l)] = receive_fg_from_neighbor(step * params.num_partitions_y + l, LEFT);

    for(uint k = 1; k < params.num_partitions_x-1; k++)
        fg_grid[get_index(k, 0)] = receive_fg_from_neighbor(step * params.num_partitions_x + k, BOTTOM);
}

void stepper_server::receive_fg_action_(uint t, vector_partition fg, direction to_dir)
{
    //direction is now opposite.
    fg_recv_buffs_[NUM_DIRECTIONS - to_dir - 1].store_received(t, std::move(fg));
}

void stepper_server::send_fg_to_neighbor(uint t, vector_partition fg, direction dir)
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
            send_uv_to_neighbor(step, uv_grid[get_index(params.num_partitions_x - 2, l)], BOTTOM_RIGHT);

        if (l == params.num_partitions_y - 2)
            send_uv_to_neighbor(step, uv_grid[get_index(1, l)], TOP_LEFT);

        send_uv_to_neighbor(step * params.num_partitions_y + l, uv_grid[get_index(1, l)], LEFT);
        send_uv_to_neighbor(step * params.num_partitions_y + l, uv_grid[get_index(params.num_partitions_x - 2, l)], RIGHT);
    }

    for (uint k = 1; k < params.num_partitions_x - 1; k++)
    {
        if (k == 1)
                send_uv_to_neighbor(step, uv_grid[get_index(k, 1)], BOTTOM_LEFT);

        if (k == params.num_partitions_x- 2)
                send_uv_to_neighbor(step, uv_grid[get_index(k, params.num_partitions_y - 2)], TOP_RIGHT);

        send_uv_to_neighbor(step * params.num_partitions_x + k, uv_grid[get_index(k, 1)], BOTTOM);
        send_uv_to_neighbor(step * params.num_partitions_x + k, uv_grid[get_index(k, params.num_partitions_y - 2)], TOP);
    }

    // RECEIVE
    for(uint l = 1; l < params.num_partitions_y - 1; l++)
    {
        uv_grid[get_index(0, l)] = receive_uv_from_neighbor(step * params.num_partitions_y + l, LEFT);
        uv_grid[get_index(params.num_partitions_x - 1, l)] = receive_uv_from_neighbor(step * params.num_partitions_y + l, RIGHT);
    }

    for(uint k = 1; k < params.num_partitions_x - 1; k++)
    {
        uv_grid[get_index(k, 0)] = receive_uv_from_neighbor(step * params.num_partitions_x + k, BOTTOM);
        uv_grid[get_index(k, params.num_partitions_y - 1)] = receive_uv_from_neighbor(step * params.num_partitions_x + k, TOP);
    }

    uv_grid[get_index(0, params.num_partitions_y - 1)] = receive_uv_from_neighbor(step, TOP_LEFT);
    uv_grid[get_index(params.num_partitions_x - 1, params.num_partitions_y -1)] = receive_uv_from_neighbor(step, TOP_RIGHT);
    uv_grid[get_index(0, 0)] = receive_uv_from_neighbor(step, BOTTOM_LEFT);
    uv_grid[get_index(params.num_partitions_x - 1, 0)] = receive_uv_from_neighbor(step, BOTTOM_RIGHT);
}

void stepper_server::receive_uv_action_(uint t, vector_partition uv, direction to_dir)
{
    //direction is now opposite.
    uv_recv_buffs_[NUM_DIRECTIONS - to_dir - 1].store_received(t, std::move(uv));
}

void stepper_server::send_uv_to_neighbor(uint t, vector_partition uv, direction dir)
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
        return vector_dummy;
}

template<typename T>
void stepper_server::print_grid(std::vector<grid::partition<T> > const& grid, const std::string message) const
{
    std::vector<std::vector<grid::partition_data<T> > > data;

    data.resize(params.num_partitions_x - 2);

    for (uint k = 1; k < params.num_partitions_x - 1; k++)
    {
        data[k-1].resize(params.num_partitions_y - 2);
        for (uint l = 1; l < params.num_partitions_y - 1; l++)
        {
            grid::partition_data<T> base = grid[get_index(k, l)].get_data(CENTER).get();
            data[k-1][l-1] = grid::partition_data<T>(base);
        }
    }

    boost::shared_ptr<hpx::lcos::local::promise<int> > p = boost::make_shared<hpx::lcos::local::promise<int> >();
    io::do_async_print(data, message, params.num_partitions_x - 2, params.num_partitions_y - 2, params.num_cells_per_partition_x, params.num_cells_per_partition_y, p);
}

void stepper_server::write_vtk(uint step) const
{
    std::vector<std::vector<scalar_data> > p_data;

    p_data.resize(params.num_partitions_x - 2);

    for (uint k = 1; k < params.num_partitions_x - 1; k++)
    {
        p_data[k-1].resize(params.num_partitions_y - 2);
        for (uint l = 1; l < params.num_partitions_y - 1; l++)
        {
            scalar_data base = p_grid[get_index(k, l)].get_data(CENTER).get();
            p_data[k-1][l-1] = scalar_data(base);
        }
    }

    std::vector<std::vector<vector_data> > uv_data;

    uv_data.resize(params.num_partitions_x - 2);

    for (uint k = 1; k < params.num_partitions_x - 1; k++)
    {
        uv_data[k-1].resize(params.num_partitions_y - 2);
        for (uint l = 1; l < params.num_partitions_y - 1; l++)
        {
            vector_data base = uv_grid[get_index(k, l)].get_data(CENTER).get();
            uv_data[k-1][l-1] = vector_data(base);
        }
    }

    boost::shared_ptr<hpx::lcos::local::promise<int> > p = boost::make_shared<hpx::lcos::local::promise<int> >();
    io::do_async_write(p_data, uv_data, params.dx, params.dy, step, params.i_max, params.j_max, params.num_partitions_x - 2, params.num_partitions_y - 2, params.num_cells_per_partition_x, params.num_cells_per_partition_y, p);

}

RealType stepper_server::compute_new_dt(std::pair<RealType, RealType> max_uv) const
{
    return c.tau * std::min(c.re/2. * 1./(1./(params.dx * params.dx) + 1./(params.dy * params.dy)), std::min(params.dx/max_uv.first, params.dy/max_uv.second));
}

}//namespace server
}//namespace stepper
