#ifndef NAST_HPX_GRID_SERVER_PARTITION_SERVER_HPP
#define NAST_HPX_GRID_SERVER_PARTITION_SERVER_HPP

#include "grid/partition_data.hpp"
#include "grid/send_buffer.hpp"
#include "grid/recv_buffer.hpp"
#include "grid/direction.hpp"

#include "io/config.hpp"

#include "util/cancellation_token.hpp"

#include <hpx/include/components.hpp>
#include <hpx/runtime/serialization/serialize.hpp>
#include <hpx/include/serialization.hpp>

namespace hpx { namespace serialization {

void serialize(input_archive& ar, std::bitset<6>& b, unsigned version)
{
    unsigned long u;
    ar >> u;

    b = std::move(std::bitset<6>(u));
}

void serialize(output_archive& ar, std::bitset<6> const& b, unsigned version)
{
    ar << b.to_ulong();
}

}}

namespace nast_hpx { namespace grid { namespace server {

char const* partition_basename = "/nast_hpx/partition/";
char const* residual_basename = "/nast/hpx/partition/residual";

/// component encapsulates partition_data, making it remotely available
struct HPX_COMPONENT_EXPORT partition_server
    : hpx::components::migration_support< hpx::components::component_base<partition_server> >
{
public:

    enum {
        U = 0,
        V,
        W,
        F,
        G,
        H,
        P,
        NUM_VARIABLES
    } variables;

    typedef hpx::serialization::serialize_buffer<double> buffer_type;
    typedef std::vector<hpx::shared_future<void> > future_vector;
    typedef std::vector<future_vector> future_grid;
    typedef std::vector<partition_data_2d<hpx::shared_future<void> > > future_grid_3d;
    typedef std::vector<future_grid_3d> future_grid_4d;

    partition_server() {}
    ~partition_server() {}

    partition_server(io::config const& cfg);

    void init();
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, init, init_action);

    /// return the sliced data appropriate for given direction
    hpx::future<triple<double> > do_timestep(double dt);
    HPX_DEFINE_COMPONENT_ACTION(partition_server, do_timestep, do_timestep_action);

    void set_left_boundary(buffer_type buffer, std::size_t step, std::size_t var)
    {
        recv_buffer_left_[var].set_buffer(buffer, step);
    }
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, set_left_boundary, set_left_boundary_action);

    void set_right_boundary(buffer_type buffer, std::size_t step, std::size_t var)
    {
        recv_buffer_right_[var].set_buffer(buffer, step);
    }
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, set_right_boundary, set_right_boundary_action);

    void set_bottom_boundary(buffer_type buffer, std::size_t step, std::size_t var)
    {
        recv_buffer_bottom_[var].set_buffer(buffer, step);
    }
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, set_bottom_boundary, set_bottom_boundary_action);

    void set_top_boundary(buffer_type buffer, std::size_t step, std::size_t var)
    {
        recv_buffer_top_[var].set_buffer(buffer, step);
    }
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, set_top_boundary, set_top_boundary_action);

    void set_front_boundary(buffer_type buffer, std::size_t step, std::size_t var)
    {
        recv_buffer_front_[var].set_buffer(buffer, step);
    }
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, set_front_boundary, set_front_boundary_action);

    void set_back_boundary(buffer_type buffer, std::size_t step, std::size_t var)
    {
        recv_buffer_back_[var].set_buffer(buffer, step);
    }
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, set_back_boundary, set_back_boundary_action);

    void set_back_left_boundary(buffer_type buffer, std::size_t step, std::size_t var)
    {
        recv_buffer_back_left_[var].set_buffer(buffer, step);
    }
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, set_back_left_boundary, set_back_left_boundary_action);

    void set_front_right_boundary(buffer_type buffer, std::size_t step, std::size_t var)
    {
        recv_buffer_front_right_[var].set_buffer(buffer, step);
    }
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, set_front_right_boundary, set_front_right_boundary_action);

    void set_bottom_right_boundary(buffer_type buffer, std::size_t step, std::size_t var)
    {
        recv_buffer_bottom_right_[var].set_buffer(buffer, step);
    }
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, set_bottom_right_boundary, set_bottom_right_boundary_action);

    void set_top_left_boundary(buffer_type buffer, std::size_t step, std::size_t var)
    {
        recv_buffer_top_left_[var].set_buffer(buffer, step);
    }
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, set_top_left_boundary, set_top_left_boundary_action);

    void set_back_bottom_boundary(buffer_type buffer, std::size_t step, std::size_t var)
    {
        recv_buffer_back_bottom_[var].set_buffer(buffer, step);
    }
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, set_back_bottom_boundary, set_back_bottom_boundary_action);

    void set_front_top_boundary(buffer_type buffer, std::size_t step, std::size_t var)
    {
        recv_buffer_front_top_[var].set_buffer(buffer, step);
    }
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, set_front_top_boundary, set_front_top_boundary_action);

    void cancel()
    {
        token.cancel();
    }
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, cancel, cancel_action);

    send_buffer<buffer_type, LEFT, set_right_boundary_action> send_buffer_left_;
    recv_buffer<buffer_type, LEFT> recv_buffer_left_[NUM_VARIABLES];

    send_buffer<buffer_type, RIGHT, set_left_boundary_action> send_buffer_right_;
    recv_buffer<buffer_type, RIGHT> recv_buffer_right_[NUM_VARIABLES];

    send_buffer<buffer_type, BOTTOM, set_top_boundary_action> send_buffer_bottom_;
    recv_buffer<buffer_type, BOTTOM> recv_buffer_bottom_[NUM_VARIABLES];

    send_buffer<buffer_type, TOP, set_bottom_boundary_action> send_buffer_top_;
    recv_buffer<buffer_type, TOP> recv_buffer_top_[NUM_VARIABLES];

    send_buffer<buffer_type, FRONT, set_back_boundary_action> send_buffer_front_;
    recv_buffer<buffer_type, FRONT> recv_buffer_front_[NUM_VARIABLES];

    send_buffer<buffer_type, BACK, set_front_boundary_action> send_buffer_back_;
    recv_buffer<buffer_type, BACK> recv_buffer_back_[NUM_VARIABLES];

    send_buffer<buffer_type, BACK_LEFT, set_front_right_boundary_action> send_buffer_back_left_;
    recv_buffer<buffer_type, BACK_LEFT> recv_buffer_back_left_[NUM_VARIABLES];

    send_buffer<buffer_type, FRONT_RIGHT, set_back_left_boundary_action> send_buffer_front_right_;
    recv_buffer<buffer_type, FRONT_RIGHT> recv_buffer_front_right_[NUM_VARIABLES];

    send_buffer<buffer_type, BOTTOM_RIGHT, set_top_left_boundary_action> send_buffer_bottom_right_;
    recv_buffer<buffer_type, BOTTOM_RIGHT> recv_buffer_bottom_right_[NUM_VARIABLES];

    send_buffer<buffer_type, TOP_LEFT, set_bottom_right_boundary_action> send_buffer_top_left_;
    recv_buffer<buffer_type, TOP_LEFT> recv_buffer_top_left_[NUM_VARIABLES];

    send_buffer<buffer_type, BACK_BOTTOM, set_front_top_boundary_action> send_buffer_back_bottom_;
    recv_buffer<buffer_type, BACK_BOTTOM> recv_buffer_back_bottom_[NUM_VARIABLES];

    send_buffer<buffer_type, FRONT_TOP, set_back_bottom_boundary_action> send_buffer_front_top_;
    recv_buffer<buffer_type, FRONT_TOP> recv_buffer_front_top_[NUM_VARIABLES];

protected:
    template<direction dir>
    void send_boundary(std::size_t step, std::size_t var, partition_data<hpx::shared_future<void> >& send_future);

    template<direction dir>
    void receive_boundary(std::size_t step, std::size_t var, future_grid_4d& recv_futures);

    template<direction dir>
    hpx::shared_future<void> get_dependency(std::size_t idx_block, std::size_t idy_block, std::size_t idz_block, future_grid_3d const& recv_futures,
    partition_data<hpx::shared_future<void> > const& calc_futures);

    void send_boundaries_U(partition_data<hpx::shared_future<void> >& send_futures, std::size_t step);
    void receive_boundaries_U(future_grid_4d& recv_futures, std::size_t step);
    void send_boundaries_V(partition_data<hpx::shared_future<void> >& send_futures, std::size_t step);
    void receive_boundaries_V(future_grid_4d& recv_futures, std::size_t step);
    void send_boundaries_W(partition_data<hpx::shared_future<void> >& send_futures, std::size_t step);
    void receive_boundaries_W(future_grid_4d& recv_futures, std::size_t step);
    void send_boundaries_F(partition_data<hpx::shared_future<void> >& send_futures, std::size_t step);
    void receive_boundaries_F(future_grid_4d& recv_futures, std::size_t step);
    void send_boundaries_G(partition_data<hpx::shared_future<void> >& send_futures, std::size_t step);
    void receive_boundaries_G(future_grid_4d& recv_futures, std::size_t step);
    void send_boundaries_H(partition_data<hpx::shared_future<void> >& send_futures, std::size_t step);
    void receive_boundaries_H(future_grid_4d& recv_futures, std::size_t step);
    void send_boundaries_P(partition_data<hpx::shared_future<void> >& send_futures, std::size_t step);
    void receive_boundaries_P(future_grid_4d& recv_futures, std::size_t step);

private:

    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned version)
    {
        ar & c & is_left_ & is_right_ & is_bottom_ & is_top_ & is_front_ & is_back_ & current & last & data_ & rhs_data_ & cell_type_data_
           & fluid_cells_ & obstacle_cells_ & cells_x_ & cells_y_ & idx_ & idy_
           & step_ & outcount_ & t_ & next_out_;
    }

    inline std::size_t get_id(std::size_t dir_x, std::size_t dir_y, std::size_t dir_z);

    partition_data<double> data_[NUM_VARIABLES];
    partition_data<double> rhs_data_;
    partition_data<std::bitset<9> > cell_type_data_;

    partition_data<std::vector<index> > fluid_cells_;
    partition_data<std::vector<index> > obstacle_cells_;

    partition_data<hpx::shared_future<void> > set_velocity_futures;
    future_grid_4d recv_futures;

    partition_data<hpx::shared_future<void> > compute_fg_futures;

    partition_data<hpx::shared_future<void> > compute_rhs_futures;

    partition_data<hpx::shared_future<double> > compute_res_futures;

    std::size_t current;
    std::size_t last;

    partition_data<hpx::shared_future<void> > set_p_futures;
    partition_data<hpx::shared_future<void> > solver_cycle_futures;

    std::vector<hpx::future<triple<double> > > local_max_uvs;

    std::vector<hpx::id_type> ids_;

    io::config c;

    std::size_t cells_x_, cells_y_, cells_z_;
    std::size_t idx_, idy_, idz_;
    std::size_t step_;
    std::size_t outcount_;

    double t_, next_out_;

    util::cancellation_token token;

    bool is_left_, is_right_, is_bottom_, is_top_, is_back_, is_front_;
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
