#ifndef GRID_SERVER_PARTITION_SERVER_HPP
#define GRID_SERVER_PARTITION_SERVER_HPP

#include <hpx/include/components.hpp>

#include "grid/partition_data.hpp"

namespace grid { namespace server {

struct HPX_COMPONENT_EXPORT partition_server
    : hpx::components::component_base<partition_server>
{
    partition_server() {}

    partition_server(partition_data<> const& data)
    : data_(data)
    {}

    partition_server(uint size_x, uint size_y, RealType initial_value)
    : data_(size_x, size_y, initial_value)
    {}

    partition_data<RealType> get_data(direction type) const
    {
        if(type == CENTER)
            return data_;
        else
            return partition_data<>(data_, type);
    }
    HPX_DEFINE_COMPONENT_ACTION(partition_server, get_data, get_data_action);


private:
    partition_data<> data_;
};

}//namespace server
}//namespace grid

HPX_REGISTER_ACTION_DECLARATION(grid::server::partition_server::get_data_action, partition_server_get_data_action);

#endif
