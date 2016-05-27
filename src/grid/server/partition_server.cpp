#include "partition_server.hpp"

#include <hpx/lcos/when_all.hpp>

typedef nast_hpx::grid::server::partition_server partition_component;
typedef hpx::components::component<partition_component> partition_server_type;

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_COMPONENT(partition_server_type, partition_component);
HPX_REGISTER_ACTION(nast_hpx::grid::server::partition_server::do_timestep_action,
    partition_server_do_timestep_action);
    
namespace nast_hpx { namespace grid { namespace server {
 
void partition_server::send_boundary()
{
    if (rank % 2 != 0)
    {
    auto f = hpx::async(hpx::util::bind(
                    boost::ref(send_buffer_left),
                    boost::ref(data_),
                    0,
                    0));
                    
    f.wait();
    }
}

void partition_server::receive_boundary()
{
    if (recv_buffer_right.valid_)
    {    
    auto f = hpx::async(hpx::util::bind(
                    boost::ref(recv_buffer_right),
                    boost::ref(data_),
                    0));

    f.wait();
    std::cout << "after got from right\n" << data_ << std::endl;
            
    }
}
    
void partition_server::do_timestep()
{
    std::cout << "ho" << std::endl;
    
    rank = hpx::get_locality_id();
    auto nl = hpx::get_num_localities_sync();
    
    std::cout << "have " << nl << " localities" << std::endl;
    
    std::vector<hpx::future<hpx::id_type > > steppers = hpx::find_all_from_basename(partition_basename, nl);
    
    ids = hpx::when_all(steppers).then(hpx::util::unwrapped2(
        [](std::vector<hpx::id_type>&& ids) -> std::vector<hpx::id_type> { return ids;})).get();
    
    std::cout << "after ids" << std::endl;
    
    if (rank % 2 != 0)
    {
        send_buffer_left.dest_ = ids[rank - 1];
    }
    else
    {
        recv_buffer_right.valid_ = true;
    }
    
    for (uint i = 0; i < data_.size_x_ ; i++)
        for (uint j = 0; j < data_.size_y_; ++j)            
            data_(i, j) = 100*rank+i+0.01*j;
            
    std::cout << data_ << std::endl;;
                 
    send_boundary();
    receive_boundary();
}   
    
}
}
}