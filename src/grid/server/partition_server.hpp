#ifndef NAST_HPX_GRID_SERVER_PARTITION_SERVER_HPP
#define NAST_HPX_GRID_SERVER_PARTITION_SERVER_HPP

#include <hpx/include/components.hpp>
#include <hpx/runtime/serialization/serialize.hpp>

#include "grid/partition_data.hpp"
#include "grid/send_buffer.hpp"
#include "grid/recv_buffer.hpp"

namespace nast_hpx { namespace grid { namespace server {

char const* partition_basename = "/nast_hpx/partition/";

    
/// component encapsulates partition_data, making it remotely available
struct HPX_COMPONENT_EXPORT partition_server
    : public hpx::components::component_base<partition_server>
{
public:
    typedef hpx::serialization::serialize_buffer<Real> buffer_type;

    partition_server() {}

    /// construct partition from given data.
    /// note: this does not copy, data_ will only point to the given data,
    /// since partition_data is essentially a shared_array
    partition_server(partition_data<Real> const& data)
    : data_(data)
    {}

    partition_server(uint size_x, uint size_y, Real initial_value)
    : data_(size_x, size_y, initial_value)
    {}

    /// return the sliced data appropriate for given direction
    void do_timestep();
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, do_timestep, do_timestep_action);
    
    void set_right_boundary(buffer_type buffer, std::size_t step, std::size_t var)
    {
        recv_buffer_right.set_buffer(buffer, step);
    }
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, set_right_boundary, set_right_boundary_action);
    
    send_buffer<buffer_type, LEFT, set_right_boundary_action> send_buffer_left;
    recv_buffer<buffer_type, RIGHT> recv_buffer_right;
        
protected:
    void send_boundary();
    void receive_boundary();

private:
    partition_data<Real> data_;
    std::vector<hpx::id_type> ids;
    std::size_t rank;
};

}//namespace server
}//namespace grid
}

// boilerplate needed

HPX_REGISTER_ACTION_DECLARATION(nast_hpx::grid::server::partition_server::do_timestep_action,
                                    partition_server_do_timestep_action);

#endif
