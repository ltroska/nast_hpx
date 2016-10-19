#ifndef NAST_HPX_GRID_SERVER_PARTITION_SERVER_HPP
#define NAST_HPX_GRID_SERVER_PARTITION_SERVER_HPP

#include <hpx/include/components.hpp>
#include <hpx/runtime/serialization/serialize.hpp>
#include <hpx/include/serialization.hpp>

#include "grid/partition_data.hpp"
#include "grid/send_buffer.hpp"
#include "grid/recv_buffer.hpp"

#include "io/config.hpp"

#include "util/cancellation_token.hpp"

namespace hpx { namespace serialization {

void serialize(input_archive& ar, std::bitset<7>& b, unsigned version)
{
    unsigned long u;
    ar >> u;

    b = std::move(std::bitset<7>(u));
}

void serialize(output_archive& ar, std::bitset<7> const& b, unsigned version)
{
    ar << b.to_ulong();
}



}}

namespace nast_hpx { namespace grid { namespace server {

char const* partition_basename = "/nast_hpx/partition/";
char const* residual_basename = "/nast/hpx/partition/residual";

/// component encapsulates partition_data, making it remotely available
struct HPX_COMPONENT_EXPORT partition_server
    : hpx::components::component_base<partition_server>
{
public:

    typedef std::vector<pair<Real> > particle_grid;

    static const std::size_t U = 0;
    static const std::size_t V = 1;
    static const std::size_t F = 2;
    static const std::size_t G = 3;
    static const std::size_t P = 4;
    static const std::size_t RHS = 5;
    static const std::size_t NUM_VARIABLES = 6;

    partition_server() {}

    partition_server(io::config const& cfg);

    void init();
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, init, init_action);

    pair<Real> do_timestep(Real dt);
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, do_timestep, do_timestep_action);

private:

    template<direction dir>
    hpx::shared_future<void> get_dependency(std::size_t idx_block, std::size_t idy_block, partition_data<hpx::shared_future<void> > calc_futures);

    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned version)
    {
        ar & c & current & last & data_ & cell_type_data_
           & fluid_cells_ & obstacle_cells_ & cells_x_ & cells_y_
           & step_ & outcount_ & t_ & next_out_;
    }

    partition_data<Real> data_[NUM_VARIABLES];
    partition_data<std::bitset<7> > cell_type_data_;
    partition_data<std::bitset<5> > empty_marker_data_;
    particle_grid particles;

    partition_data<std::vector<index> > fluid_cells_;
    partition_data<std::vector<index> > obstacle_cells_;

    partition_data<hpx::shared_future<void> > set_velocity_futures;
    partition_data<hpx::shared_future<void> > update_particle_futures;
    partition_data<hpx::shared_future<void> > compute_fg_futures;
    partition_data<hpx::shared_future<void> > compute_rhs_futures;
    partition_data<hpx::shared_future<Real> > compute_res_futures;
    partition_data<hpx::shared_future<void> > set_p_futures;
    std::array<partition_data<hpx::shared_future<void> >, 2> solver_cycle_futures;

    std::vector<hpx::future<pair<Real> > > local_max_uvs;
    io::config c;

    std::size_t cells_x_, cells_y_;
    std::size_t step_;
    std::size_t outcount_;
    std::size_t current, last;

    Real t_, next_out_;

    util::cancellation_token token;
};

}//namespace server
}//namespace grid
}

// boilerplate needed

HPX_REGISTER_ACTION_DECLARATION(nast_hpx::grid::server::partition_server::do_timestep_action,
                                    partition_server_do_timestep_action);

HPX_REGISTER_ACTION_DECLARATION(nast_hpx::grid::server::partition_server::init_action,
                                    partition_server_init_action);

#endif
