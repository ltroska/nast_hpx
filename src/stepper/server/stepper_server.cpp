#include <hpx/include/iostreams.hpp>
#include <cmath>

#include "stepper_server.hpp"
#include "io/vtk_writer.hpp"

typedef stepper::server::stepper_server stepper_component;
typedef hpx::components::component<stepper_component> stepper_server_type;

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_COMPONENT(stepper_server_type, stepper_component);

HPX_REGISTER_ACTION(stepper::server::stepper_server::do_work_action, stepper_server_do_work_action);
HPX_REGISTER_ACTION(stepper::server::stepper_server::setup_action, stepper_server_setup_action);


namespace stepper { namespace server {

uint get_neighbor_id(uint id, direction dir, uint num_localities)
{
    uint res_x, res_y;
    if (num_localities == 2)
    {
        res_x = 2;
        res_y = 1;
    }
    else
    {
        res_x = static_cast<uint>(sqrt(num_localities));
        res_y = res_x;
    }

    switch (dir)
    {
        case left:
            return ((id-1)/res_x == id/res_x && id-1 < num_localities) ? id-1 : num_localities;

        case right:
            return ((id+1)/res_x == id/res_x && id+1 < num_localities) ? id+1 : num_localities;

        case top:
            return ((id+res_x) < num_localities && (id+res_x)/res_x == id/res_x +1) ? id+res_x : num_localities;

        case bottom:
            return ((id-res_x) < num_localities && (id-res_x)/res_x == id/res_x -1) ? id-res_x : num_localities;

        case top_left:
            return ((id+res_x-1) < num_localities && (id+res_x-1)/res_x == id/res_x +1) ? id+res_x-1 : num_localities;

        case top_right:
            return ((id+res_x+1) < num_localities && (id+res_x+1)/res_x == id/res_x+1) ? id+res_x+1 : num_localities;

        case bottom_left:
            return ((id-res_x-1) < num_localities && (id-res_x-1)/res_x == id/res_x-1) ? id-res_x-1 : num_localities;

        case bottom_right:
            return ((id-res_x+1) < num_localities && (id-res_x+1)/res_x == id/res_x-1) ? id-res_x+1 : num_localities;
        default:
            return num_localities;
    }
}

//get type of stepper
stepper_type get_type(uint id, uint num_localities)
{
    uint res = static_cast<uint>(sqrt(num_localities));

    uint id_mod_res = id % res;
    uint id_div_res = id / res;

    if (id_mod_res == 0)
    {
        if (id_div_res == 0)
            return stepper_type::bottom_left_boundary;
        if (id_div_res == res-1)
            return stepper_type::top_left_boundary;
        return stepper_type::left_boundary;
    }

    if (id_mod_res == res-1)
    {
        if (id_div_res == 0)
            return stepper_type::bottom_right_boundary;
        if (id_div_res == res-1)
            return stepper_type::top_right_boundary;
        return stepper_type::right_boundary;
    }

    if (id_div_res == 0)
        return stepper_type::bottom_boundary;

    if (id_div_res == res-1)
        return stepper_type::top_boundary;

    return stepper_type::interior;
}

stepper_server::stepper_server(uint num_localities)
    : num_localities_(num_localities),
      locality_id_(hpx::get_locality_id()),
      left_(hpx::find_from_basename(stepper_basename, get_neighbor_id(hpx::get_locality_id(), left, num_localities))),
      right_(hpx::find_from_basename(stepper_basename, get_neighbor_id(hpx::get_locality_id(), right, num_localities))),
      top_(hpx::find_from_basename(stepper_basename, get_neighbor_id(hpx::get_locality_id(), top, num_localities))),
      bottom_(hpx::find_from_basename(stepper_basename, get_neighbor_id(hpx::get_locality_id(), bottom, num_localities))),
      top_left_(hpx::find_from_basename(stepper_basename, get_neighbor_id(hpx::get_locality_id(), top_left, num_localities))),
      top_right_(hpx::find_from_basename(stepper_basename, get_neighbor_id(hpx::get_locality_id(), top_right, num_localities))),
      bottom_left_(hpx::find_from_basename(stepper_basename, get_neighbor_id(hpx::get_locality_id(), bottom_left, num_localities))),
      bottom_right_(hpx::find_from_basename(stepper_basename, get_neighbor_id(hpx::get_locality_id(), bottom_right, num_localities))),
      type_(get_type(hpx::get_locality_id(), num_localities))
{
    hpx::cout << "on " << hpx::get_locality_id() << " left " << get_neighbor_id(hpx::get_locality_id(), left, num_localities) << " right " << get_neighbor_id(hpx::get_locality_id(), right, num_localities)
    << " top " << get_neighbor_id(hpx::get_locality_id(), top, num_localities) << " bottom " << get_neighbor_id(hpx::get_locality_id(), bottom, num_localities)
    << " top left "<<get_neighbor_id(hpx::get_locality_id(), top_left, num_localities) << " top right " << get_neighbor_id(hpx::get_locality_id(), top_right, num_localities)
    << " bottom left " << get_neighbor_id(hpx::get_locality_id(), bottom_left, num_localities) << " bottom right " << get_neighbor_id(hpx::get_locality_id(), bottom_right, num_localities)<< hpx::endl << hpx::flush;
}

uint stepper_server::setup(uint i_max, uint j_max, RealType x_length, RealType y_length, uint num_partitions_x, uint num_partitions_y)
{
    //padding for neighbor locality values
    num_local_partitions_x_ = num_partitions_x + 2;
    num_local_partitions_y_ = num_partitions_y + 2;

    x_length_ = x_length;
    y_length_ = y_length;

    dx_ = 1./i_max;
    dy_ = 1./j_max;

    //how many steppers for each dimension (i.e. if num_localities_ == 4, we need 2x2 steppers).
    if (num_localities_ == 2)
    {
        res_x_ = 2;
        res_y_ = 1;
    }
    else
    {
        res_x_ = static_cast<uint>(sqrt(num_localities_));
        res_y_ = res_x_;
    }

    num_cells_x_ = (i_max / num_localities_) / num_partitions_x;
    num_cells_y_ = (j_max / num_localities_) / num_partitions_x;

    U.resize(num_local_partitions_x_);
    for (uint i = 0; i < num_local_partitions_x_; ++i)
        U[i].resize(num_local_partitions_y_);

    uint num_local_cells_x, num_local_cells_y;

    for (uint j = 0; j < num_local_partitions_y_; ++j)
    {
        for (uint i = 0; i < num_local_partitions_x_; ++i)
        {
            num_local_cells_x = num_cells_x_;
            num_local_cells_y = num_cells_y_;
            if (i == 0 || i == num_local_partitions_x_ - 1)
                num_local_cells_x = 1;
            if (j == 0 || j == num_local_partitions_y_ - 1)
                num_local_cells_y = 1;
            U[i][j] = grid::partition(hpx::find_here(), num_local_cells_x, num_local_cells_y,
                                    ( i == 0 || j == 0 || i == num_local_partitions_x_ -1 || j == num_local_partitions_y_ -1 ) ? hpx::get_locality_id() +4 : hpx::get_locality_id() + j*0.1 + i*0.01);
        }
    }

   if (get_neighbor_id(hpx::get_locality_id(), left, num_localities_) < num_localities_)
    {
            hpx::cout << "sending left on " << hpx::get_locality_id() << hpx::endl << hpx::flush;

        for (uint j = 1; j < num_local_partitions_y_-1; j++)
            send_to_neighbor(j, U[1][j], left);
    }

    if (get_neighbor_id(hpx::get_locality_id(), right, num_localities_) < num_localities_)
    {
        hpx::cout << "sending right on " << hpx::get_locality_id() << hpx::endl << hpx::flush;
        for (uint j = 1; j < num_local_partitions_y_-1; j++)
            send_to_neighbor(j, U[num_local_partitions_x_-2][j], right);
    }

    if (get_neighbor_id(hpx::get_locality_id(), bottom, num_localities_) < num_localities_)
    {
            hpx::cout << "sending bottom on " << hpx::get_locality_id() << hpx::endl << hpx::flush;

        for (uint i = 1; i < num_local_partitions_x_-1; i++)
            send_to_neighbor(i, U[i][1], bottom);
    }

    if (get_neighbor_id(hpx::get_locality_id(), top, num_localities_) < num_localities_)
    {
            hpx::cout << "sending top on " << hpx::get_locality_id() << hpx::endl << hpx::flush;

        for (uint i = 1; i < num_local_partitions_x_-1; i++)
            send_to_neighbor(i, U[i][num_local_partitions_y_-2], top);
    }

    if (get_neighbor_id(hpx::get_locality_id(), top_left, num_localities_) < num_localities_)
    {
            hpx::cout << "sending topleft on " << hpx::get_locality_id() << hpx::endl << hpx::flush;

        send_to_neighbor(0, U[1][num_local_partitions_y_-2], top_left);
    }

    if (get_neighbor_id(hpx::get_locality_id(), top_right, num_localities_) < num_localities_)
    {
            hpx::cout << "sending topright on " << hpx::get_locality_id() << hpx::endl << hpx::flush;

        send_to_neighbor(0, U[num_local_partitions_x_-2][num_local_partitions_y_-2], top_right);
    }

    if (get_neighbor_id(hpx::get_locality_id(), bottom_left, num_localities_) < num_localities_)
    {
            hpx::cout << "sending bottomleft on " << hpx::get_locality_id() << hpx::endl << hpx::flush;

        send_to_neighbor(0, U[1][1], bottom_left);
    }

    if (get_neighbor_id(hpx::get_locality_id(), bottom_right, num_localities_) < num_localities_)
    {
            hpx::cout << "sending bottomright on " << hpx::get_locality_id() << hpx::endl << hpx::flush;

        send_to_neighbor(0, U[1][num_local_partitions_y_-2], bottom_right);
    }

    hpx::cout << "stepper on " << locality_id_ << " with " << num_cells_x_ << "x" << num_cells_y_ << " cells, " << num_local_partitions_x_
    << "x"<< num_local_partitions_y_ << " partitions, dx=" << dx_ << " dy=" << dy_ << hpx::endl << hpx::flush;
}

uint stepper_server::do_work(uint num_local_partitions_x, uint num_local_partitions_y, uint num_cells_x, uint num_cells_y, RealType delta_x, RealType delta_y)
{
    write_vtk_files();
    return 0;
}

void stepper_server::set_boundary_values_u_v()
{
    for (auto col : U)
    {
        for (auto part : col)
        {
            grid::partition_data<grid::scalar_cell> pdata = part.get_p_data(grid::center_partition).get();

            for(uint i = 0; i < pdata.size_x(); ++i)
                for(uint j = 0; j < pdata.size_y(); ++j)
                {
                    grid::scalar_cell& c = pdata.get_cell_ref(i,j);
                    c.c = hpx::get_locality_id();
                }
        }
    }
}

void stepper_server::write_vtk_files()
{
    std::vector<std::vector<grid::partition_data<grid::scalar_cell> > > space_part;
    space_part.resize(num_local_partitions_x_);

    for (uint i = 0; i < num_local_partitions_x_; ++i)
    {
        space_part[i].resize(num_local_partitions_y_);
        for (uint j = 0; j < num_local_partitions_y_; ++j)
        {
            grid::partition_data<grid::scalar_cell> base = U[i][j].get_p_data(grid::partition_type::center_partition).get();
            space_part[i][j] = grid::partition_data<grid::scalar_cell>(base.size_x(), base.size_y());
            for (uint k = 0; k < base.size(); ++k) {
                space_part[i][j][k] = base[k];
            }
        }
    }

    io::do_async_write(space_part, num_local_partitions_x_, num_local_partitions_y_, num_cells_x_, num_cells_y_);
}

}//namespace server
}//namespace stepper
