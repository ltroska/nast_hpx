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
    strategy->set_velocity_on_boundary(uv_grid);
    print_grid(uv_grid, "uv_data");

    std::cout << "stepper on " << hpx::get_locality_id() << " with " <<  num_partitions_x-2 << "x"<< num_partitions_y-2 << " partitions, "
             << params.num_cells_per_partition_x << "x" <<  params.num_cells_per_partition_y << " cells each, " << "dx=" << params.dx << " dy=" << params.dx
             << std::endl;
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

    num_partitions_x = c.i_res + 2;
    num_partitions_y = c.j_res + 2;

    params.re = c.re;
    params.alpha = c.alpha;
    params.omega = c.omega;
    params.dx = c.x_length / c.i_max;
    params.dy = c.y_length / c.j_max;

    params.i_max = c.i_max;
    params.j_max = c.j_max;
    params.num_partitions_x = num_partitions_x;
    params.num_partitions_y = num_partitions_y;
    params.num_cells_per_partition_x = ((c.i_max + 2) / num_localities_x) / c.i_res;
    params.num_cells_per_partition_y = ((c.j_max + 2) / num_localities_y) / c.j_res;
}

void stepper_server::initialize_grids()
{
    index_grid.resize(num_partitions_x * num_partitions_y);
    uv_grid.resize(num_partitions_x * num_partitions_y);

    for (uint l = 0; l < num_partitions_y; l++)
        for (uint k = 0; k < num_partitions_x; k++)
        {
            index_grid[get_index(k, l)] =
                std::pair<RealType, RealType>
                    (
                    (hpx::get_locality_id() % num_localities_x) * (num_partitions_x - 2) * params.num_cells_per_partition_x + (k - 1) * params.num_cells_per_partition_x,
                    (hpx::get_locality_id() / num_localities_x) * (num_partitions_y - 2) * params.num_cells_per_partition_y + (l - 1) * params.num_cells_per_partition_y
                    );

            uv_grid[get_index(k, l)] = vector_partition(hpx::find_here(), params.num_cells_per_partition_x, params.num_cells_per_partition_y);
        }
}

uint stepper_server::get_index(uint k, uint l) const
{
    return l * num_partitions_x + k;
}

template<typename T>
void stepper_server::print_grid(std::vector<grid::partition<T> > const& grid, const std::string message) const
{
    std::vector<std::vector<grid::partition_data<T> > > data;

    data.resize(num_partitions_x - 2);

    for (uint k = 1; k < num_partitions_x - 1; k++)
    {
        data[k-1].resize(num_partitions_y - 2);
        for (uint l = 1; l < num_partitions_y - 1; l++)
        {
            grid::partition_data<T> base = grid[get_index(k, l)].get_data(CENTER).get();
            data[k-1][l-1] = grid::partition_data<T>(base);
        }
    }

    boost::shared_ptr<hpx::lcos::local::promise<int> > p = boost::make_shared<hpx::lcos::local::promise<int> >();
    io::do_async_print(data, message, num_partitions_x - 2, num_partitions_y - 2, params.num_cells_per_partition_x, params.num_cells_per_partition_y, p);
}

}//namespace server
}//namespace stepper
