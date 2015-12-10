#ifndef GRID_SERVER_PARTITION_SERVER_HPP
#define GRID_SERVER_PARTITION_SERVER_HPP

#include <hpx/include/components.hpp>

#include "grid/partition_data.hpp"
#include "internal/types.hpp"

namespace grid { namespace server {

struct HPX_COMPONENT_EXPORT partition_server
    : hpx::components::component_base<partition_server>
{
    partition_server() {}

    partition_server(partition_data<scalar_cell> const& p_data)
    : p_data_(p_data)
    {}

    partition_server(uint size_x, uint size_y, RealType initial_value)
    : p_data_(size_x, size_y, initial_value)
    {}

    partition_data<scalar_cell> get_p_data(partition_type type) const
    {
        if(type == center_partition)
            return p_data_;
        else
            return partition_data<scalar_cell>(p_data_, type);
    }

    HPX_DEFINE_COMPONENT_ACTION(partition_server, get_p_data, get_p_data_action);

    partition_data<vector_cell> get_uv_data(partition_type type) const
    {
        if(type == center_partition)
            return uv_data_;
        else
            return partition_data<vector_cell>(uv_data_, type);
    }

    HPX_DEFINE_COMPONENT_ACTION(partition_server, get_uv_data, get_uv_data_action);

    partition_data<vector_cell> get_fg_data(partition_type type) const
    {
        if(type == center_partition)
            return fg_data_;
        else
            return partition_data<vector_cell>(fg_data_, type);
    }

    HPX_DEFINE_COMPONENT_ACTION(partition_server, get_fg_data, get_fg_data_action);

private:
    partition_data<scalar_cell> p_data_;
    partition_data<vector_cell> uv_data_;
    partition_data<vector_cell> fg_data_;
};

}//namespace server
}//namespace grid

HPX_REGISTER_ACTION_DECLARATION(grid::server::partition_server::get_p_data_action, partition_server_get_p_data_action);
HPX_REGISTER_ACTION_DECLARATION(grid::server::partition_server::get_uv_data_action, partition_server_get_uv_data_action);
HPX_REGISTER_ACTION_DECLARATION(grid::server::partition_server::get_fg_data_action, partition_server_get_fg_data_action);
#endif
