#include <hpx/include/iostreams.hpp>
#include <cmath>

#include "stepper_server.hpp"
#include "io/vtk_writer.hpp"

typedef stepper::server::stepper_server stepper_component;
typedef hpx::components::component<stepper_component> stepper_server_type;

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_COMPONENT(stepper_server_type, stepper_component);

HPX_REGISTER_ACTION(stepper::server::stepper_server::do_work_action, stepper_server_do_work_action);


namespace stepper { namespace server {

enum direction
{
    left, right, top, bottom, top_left, top_right, bottom_left, bottom_right
};

enum type
{
    left_boundary, right_boundary, top_boundary, bottom_boundary, top_left_boundary, top_right_boundary, bottom_left_boundary, bottom_right_boundary, interior
};

uint get_neighbor_id(uint id, direction dir, uint num_localities)
{
    uint res = static_cast<uint>(sqrt(num_localities));

    switch (dir)
    {
        case left:
            return ((id-1)/res == id/res && id-1 < num_localities) ? id-1 : num_localities;

        case right:
            return ((id+1)/res == id/res && id+1 < num_localities) ? id+1 : num_localities;

        case top:
            return ((id+res) < num_localities && (id+res)/res == id/res +1) ? id+res : num_localities;

        case bottom:
            return ((id-res) < num_localities && (id-res)/res == id/res -1) ? id-res : num_localities;

        case top_left:
            return ((id+res-1) < num_localities && (id+res-1)/res == id/res +1) ? id+res-1 : num_localities;

        case top_right:
            return ((id+res+1) < num_localities && (id+res+1)/res == id/res+1) ? id+res+1 : num_localities;

        case bottom_left:
            return ((id-res-1) < num_localities && (id-res-1)/res == id/res-1) ? id-res-1 : num_localities;

        case bottom_right:
            return ((id-res+1) < num_localities && (id-res+1)/res == id/res-1) ? id-res+1 : num_localities;
        default:
            return num_localities;
    }
}

//get type of stepper
type get_type(uint id, uint num_localities)
{
    uint res = static_cast<uint>(sqrt(num_localities));

    uint id_mod_res = id % res;
    uint id_div_res = id / res;

    if (id_mod_res == 0)
    {
        if (id_div_res == 0)
            return type::bottom_left_boundary;
        if (id_div_res == res-1)
            return type::top_left_boundary;
        return type::left_boundary;
    }

    if (id_mod_res == res-1)
    {
        if (id_div_res == 0)
            return type::bottom_right_boundary;
        if (id_div_res == res-1)
            return type::top_right_boundary;
        return type::right_boundary;
    }

    if (id_div_res == 0)
        return type::bottom_boundary;

    if (id_div_res == res-1)
        return type::top_boundary;

    return type::interior;
}

stepper_server::stepper_server(uint num_localities)
    : num_localities_(num_localities),
      left_(hpx::find_from_basename(stepper_basename,0))
{
   hpx::cout << "new stepper on locality " << hpx::get_locality_id() << hpx::endl << hpx::flush;
}

uint stepper_server::do_work(uint num_local_partitions_x, uint num_local_partitions_y, uint num_cells_x, uint num_cells_y, RealType delta_x, RealType delta_y)
{

    num_local_partitions_x_ = num_local_partitions_x;
    num_local_partitions_y_ = num_local_partitions_y;
    num_cells_x_ = num_cells_x;
    num_cells_y_ = num_cells_y;
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
  // write_vtk_files();

    if (hpx::get_locality_id() == 1 || hpx::get_locality_id() == 3)
    {
        grid::partition p = U[0][0];
        grid::partition_data pdata = p.get_data(grid::partition_type::center_partition).get();
    send_left(0, U[0][0]);

    }

    if (hpx::get_locality_id() == 0 || hpx::get_locality_id() == 2)
    {
        grid::partition p = receive_right(0);
        grid::partition_data pdata = p.get_data(grid::partition_type::top_left_partition).get();

    }
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
    space_part.resize(num_local_partitions_x_);

    for (uint i = 0; i < num_local_partitions_x_; ++i)
    {
        space_part[i].resize(num_local_partitions_y_);
        for (uint j = 0; j < num_local_partitions_y_; ++j)
        {
            grid::partition_data base = U[i][j].get_data(grid::partition_type::center_partition).get();
            space_part[i][j] = grid::partition_data(base.size_x(), base.size_y());
            for (uint k = 0; k < base.size(); ++k) {
                space_part[i][j][k] = base[k];
            }
        }
    }

    io::do_async_write(space_part, num_local_partitions_x_, num_local_partitions_y_, num_cells_x_, num_cells_y_);
}

}//namespace server
}//namespace stepper
