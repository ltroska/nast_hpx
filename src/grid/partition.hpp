#ifndef GRID_PARTITION_HPP
#define GRID_PARTITION_HPP

#include "server/partition_server.hpp"

namespace grid {

/// This class is a client for the partition_server component, simplifying
/// access.
template <typename T=RealType>
struct partition
    : hpx::components::client_base<partition<T>, server::partition_server<T> >
{
    typedef hpx::components::client_base<partition<T>,
                                            server::partition_server<T> >
        base_type;
    
    typedef server::partition_server<T> server_type;

    partition() {}

    partition(hpx::id_type where, uint size_x, uint size_y,
                RealType initial_value)
        : base_type(hpx::new_<server::partition_server<T> >(where, size_x,
            size_y, initial_value))
    {}

    partition(hpx::id_type where, uint size_x, uint size_y)
        : base_type(hpx::new_<server::partition_server<T> >(where, size_x,
            size_y, 0))
    {
    }

    // Create a new component on the locality co-located to the id 'where'. The
    // new instance will be initialized from the given partition_data.
    partition(hpx::id_type where, partition_data<T> const& data)
      : base_type(hpx::new_<server::partition_server<T> >(hpx::colocated(where),
            data))
    {
    }

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
    {
    }

    // Get the sliced data from the component.
    hpx::future<partition_data<T> > get_data(direction type) const
    {
        typename server_type::get_data_action act;
        return hpx::async(act, this->get_id(), type);
    }
};
}//namespace grid
#endif
