#include <hpx/include/iostreams.hpp>

#include <cmath>

#include "stepper_server.hpp"
#include "computation/custom_grain_size.hpp"

#include "io/vtk_writer.hpp"

typedef stepper::server::stepper_server stepper_component;
typedef hpx::components::component<stepper_component> stepper_server_type;

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_COMPONENT(stepper_server_type, stepper_component);

namespace stepper { namespace server {

stepper_server::stepper_server(uint nl)
    : num_localities(nl)
{}

stepper_server::stepper_server(uint nl, io::config const& cfg)
    : num_localities(nl), c(cfg)
{
    initialize_parameters();
    initialize_grids();

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
}

void stepper_server::do_work()
{
    std::pair<RealType, RealType> max_uv(2, 0);

    RealType t = 0, dt = 0;
    for (uint step = 0; t < c.t_end; step++)
    {
        dt = compute_new_dt(max_uv);
        t += dt;

        std::cout << "step: " << step << " t: " << t << " dt: " << dt;
        max_uv = do_timestep(step, dt).get();
        std:: cout << " max_uv: " << max_uv.first << " " << max_uv.second << std::endl;

       // print_grid(p_grid);
       if (step % c.output_skip_size == 0)
            write_vtk(step);
    }
}

hpx::future<std::pair<RealType, RealType> > stepper_server::do_timestep(uint step, RealType dt)
{
    strategy->set_velocity_on_boundary(uv_grid);

    strategy->compute_fg(fg_grid, uv_grid, dt);

    strategy->compute_rhs(rhs_grid, fg_grid, dt);

    RealType residual;
    uint iter;
    for (iter = 0; iter < c.iter_max; iter++)
    {
        strategy->set_pressure_on_boundary(p_grid);

        residual = strategy->sor_cycle(p_grid, rhs_grid).get();


        if (residual < c.eps_sq)
            break;
    }

    std::cout << " iterations: " << iter << " residual: " << residual;

    return strategy->update_velocities(uv_grid, fg_grid, p_grid, dt);
}

uint stepper_server::get_index(uint k, uint l) const
{
    return l * params.num_partitions_x + k;
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

void stepper_server::write_vtk(uint step)
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

RealType stepper_server::compute_new_dt(std::pair<RealType, RealType> max_uv)
{
    return c.tau * std::min(c.re/2. * 1./(1./(params.dx * params.dx) + 1./(params.dy * params.dy)), std::min(params.dx/max_uv.first, params.dy/max_uv.second));
}

}//namespace server
}//namespace stepper
