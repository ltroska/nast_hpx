#ifndef GRID_SERVER_PARTITION_SERVER_HPP
#define GRID_SERVER_PARTITION_SERVER_HPP

#include <hpx/include/iostreams.hpp>
#include <hpx/include/components.hpp>

#include "grid/partition_data.hpp"
#include "internal/types.hpp"

namespace grid { namespace server {

struct HPX_COMPONENT_EXPORT partition_server
    : hpx::components::component_base<partition_server>
{
    partition_server() {}

    partition_server(partition_data const& data)
    : data_(data)
    {}

    partition_server(uint size_x, uint size_y, RealType initial_value)
    : data_(size_x, size_y, initial_value)
    {}

    partition_data get_data(partition_type type) const
    {
        if(type == center_partition)
            return data_;
        else
            return partition_data(data_, type);
    }

    uint set_boundary(boundary_type type, RealType constant_velocity)
    {
        uint size_x = data_.size_x();
        uint size_y = data_.size_y();


        switch (type)
        {
            case left_boundary:
            {
                for (uint j = 0; j != size_y; ++j)
                {
                    cell& current_cell = data_.get_cell_ref(0, j);
                    cell neighbor_cell = data_.get_cell(1,j);
                    /*
                    *@todo this is intentionally wrong
                    */
                    current_cell.p = (constant_velocity == 0) ? neighbor_cell.p : constant_velocity;
                    current_cell.u = 0;
                    current_cell.v = -neighbor_cell.v;

                    current_cell.F = current_cell.u;
                }

                break;
            }

            case right_boundary:
            {
                for (uint j = 0; j != size_y; ++j)
                {
                    cell& current_cell = data_.get_cell_ref(size_x-1, j);
                    cell& neighbor_cell = data_.get_cell_ref(size_x-2,j);

                    neighbor_cell.u = 0;

                    current_cell.p = (constant_velocity == 0) ? neighbor_cell.p : constant_velocity;
                    current_cell.v = -neighbor_cell.v;

                    neighbor_cell.F = neighbor_cell.u;
                }

                break;
            }

            case bottom_boundary:
            {
                for (uint i = 0; i != size_x; ++i)
                {
                    cell& current_cell = data_.get_cell_ref(i, 0);
                    cell neighbor_cell = data_.get_cell(i,1);

                    current_cell.p = (constant_velocity == 0) ? neighbor_cell.p : constant_velocity;
                    current_cell.u = -neighbor_cell.u;
                    current_cell.v = 0;

                    current_cell.G = current_cell.v;
                }

                break;
            }

            case top_boundary:
            {
                for (uint i = 0; i != size_x; ++i)
                {
                    cell& current_cell = data_.get_cell_ref(i, size_y-1);
                    cell& neighbor_cell = data_.get_cell_ref(i,size_y-2);

                    neighbor_cell.v = 0;

                    current_cell.p = (constant_velocity == 0) ? neighbor_cell.p : constant_velocity;
                    current_cell.u = -neighbor_cell.u;

                    neighbor_cell.G = neighbor_cell.v;
                }

                break;
            }

            default:
                break;
        }

        return 1;
    }

    HPX_DEFINE_COMPONENT_ACTION(partition_server, get_data, get_data_action);
    HPX_DEFINE_COMPONENT_ACTION(partition_server, set_boundary, set_boundary_action);

private:
    partition_data data_;
};

}//namespace server
}//namespace grid

HPX_REGISTER_ACTION_DECLARATION(grid::server::partition_server::get_data_action, partition_server_get_data_action);
HPX_REGISTER_ACTION_DECLARATION(grid::server::partition_server::set_boundary_action, partition_server_set_boundary_action);
#endif
