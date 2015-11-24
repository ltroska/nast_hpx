#ifndef GRID_PARTITION_HPP
#define GRID_PARTITION_HPP

#include <hpx/hpx.hpp>

#include "server/partition_server.hpp"
#include "internal/types.hpp"

namespace grid {

struct partition
    : hpx::components::client_base<partition, server::partition_server>
{
    typedef hpx::components::client_base<partition, server::partition_server> base_type;

    partition() {}

    partition(hpx::id_type where, uint size_x, uint size_y, RealType initial_value)
        : base_type(hpx::new_<server::partition_server>(where, size_x, size_y, initial_value))
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

    hpx::future<partition_data> get_data(partition_type type) const
    {
        server::partition_server::get_data_action act;
        return hpx::async(act, get_id(), type);
    }

};

}//namespace grid

#endif
