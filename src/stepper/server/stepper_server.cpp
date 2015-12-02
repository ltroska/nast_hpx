#include <hpx/include/iostreams.hpp>

#include "stepper_server.hpp"
#include "io/vtk_writer.hpp"

typedef stepper::server::stepper_server stepper_component;
typedef hpx::components::component<stepper_component> stepper_server_type;

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_COMPONENT(stepper_server_type, stepper_component);

HPX_REGISTER_ACTION(stepper::server::stepper_server::do_work_action, stepper_server_do_work_action);


namespace stepper { namespace server {
//calculate locality to distribute partitions evenly
inline uint locidx(uint i, uint np, uint nl)
{
    return i / (np/nl);
}

stepper_server::stepper_server()
{
   hpx::cout << "new stepper on locality " << hpx::get_locality_id() << hpx::endl << hpx::flush;
}

stepper_server::stepper_server(uint num_partitions_x, uint num_partitions_y, uint cells_x, uint cells_y, RealType delta_x, RealType delta_y)
    : num_local_partitions_x(num_partitions_x), num_local_partitions_y(num_partitions_y), num_cells_x(cells_x), num_cells_y(cells_y), dx(delta_x), dy(delta_y)
{
   hpx::cout << "new stepper on locality " << hpx::get_locality_id() << hpx::endl << hpx::flush;
}

uint stepper_server::do_work()
{
    U.resize(num_local_partitions_x);

    for (uint i = 0; i < num_local_partitions_x; ++i)
        U[i].resize(num_local_partitions_y);

    for (uint j = 0; j < num_local_partitions_y; ++j)
    {
        for (uint i = 0; i < num_local_partitions_x; ++i)
        {
            U[i][j] = grid::partition(hpx::find_here(), num_cells_x, num_cells_y, hpx::get_locality_id() + j*num_local_partitions_x + i);
        }
    }

    set_boundary_values_u_v();
    write_vtk_files();

    // hpx::cout << hpx::find_here() << hpx::flush()
    return 0;
}

void stepper_server::set_boundary_values_u_v()
{
    for (auto col : U)
    {
        for (auto part : col)
        {
            grid::partition_data pdata = part.get_data(grid::center_partition).get();

            for(uint i = 0; i < pdata.size_x(); ++i)
                for(uint j = 0; j < pdata.size_y(); ++j)
                {
                    grid::cell& c = pdata.get_cell_ref(i,j);
                    c.p = hpx::get_locality_id();
                }
        }
    }
}

void stepper_server::write_vtk_files()
{
    std::vector<std::vector<grid::partition_data> > space_part;
    space_part.resize(num_local_partitions_y);

    for (uint i = 0; i < num_local_partitions_x; ++i)
        for (uint j = 0; j < num_local_partitions_y; ++j)
            space_part[i].push_back(U[i][j].get_data(grid::center_partition).get());


    io::async_write(space_part, num_local_partitions_x, num_local_partitions_y, num_cells_x, num_cells_y);
}

}//namespace server
}//namespace stepper
