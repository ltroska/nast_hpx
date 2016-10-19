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

    partition(hpx::id_type where, io::config const& cfg)
        : base_type(hpx::new_<server::partition_server>(where, cfg))
    {}

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

    partition(hpx::shared_future<hpx::id_type> const& id)
      : base_type(id)
    {}

    partition(partition const& other)
      : base_type(other)
    {}

    hpx::future<void> init()
    {
        return hpx::async<server::partition_server::init_action>(this->get_id());
    }

    void init_sync()
    {
        init().wait();
    }

    hpx::future<pair<Real> > do_timestep(Real dt)
    {
        typename server::partition_server::do_timestep_action act;
        return hpx::async(act, get_id(), dt);
    }

    hpx::future<pair<Real> > do_timestep_fut(hpx::future<Real> dt)
    {
        typename server::partition_server::do_timestep_action act;
        return hpx::dataflow(act, get_id(), dt.get());
    }
};

}//namespace grid
}

#endif
