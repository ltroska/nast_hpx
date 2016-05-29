#ifndef NAST_HPX_GRID_PARTITION_HPP
#define NAST_HPX_GRID_PARTITION_HPP

#include "server/partition_server.hpp"

namespace nast_hpx { namespace grid {

/// This class is a client for the partition_server component, simplifying
/// access.
struct partition
    : hpx::components::client_base<partition, server::partition_server>
{
    typedef hpx::components::client_base<partition, server::partition_server>
        base_type;
    
    partition() {}

    partition(hpx::id_type where, io::config& cfg, std::size_t idx, std::size_t idy)
        : base_type(hpx::new_<server::partition_server>(where, cfg, idx, idy))
    {
        std::cout << "registering " << idy * cfg.num_localities_x + idx << std::endl;
        hpx::register_with_basename(server::partition_basename, get_id(),
                                        idy * cfg.num_localities_x + idx);    
    }

    // Create a new component on the locality co-located to the id 'where'. The
    // new instance will be initialized from the given partition_data.
    partition(hpx::id_type where)
      : base_type(hpx::new_<server::partition_server>(hpx::colocated(where)))
    {}

    // Attach a future representing a (possibly remote) partition.
    partition(hpx::future<hpx::id_type> && id)
      : base_type(std::move(id))
    {}

    // Unwrap a future<partition> (a partition already holds a future to the
    // id of the referenced object, thus unwrapping accesses this inner future).
    partition(hpx::future<partition> && c)
      : base_type(std::move(c))
    {}

    partition(partition const& other)
      : base_type(other)
    {}
    
    hpx::future<void> init()
    {
        typename server::partition_server::init_action act;
        return hpx::async(act, get_id());  
    }    
    
    void init_sync()
    {
        typename server::partition_server::init_action act;
        hpx::apply(act, get_id());  
    }

    hpx::future<void> do_timestep(Real dt) 
    {
        typename server::partition_server::do_timestep_action act;
        return hpx::async(act, get_id(), dt);
    }
};

}//namespace grid
}

#endif
