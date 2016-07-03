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

    static const std::size_t U = 0;
    static const std::size_t V = 1;
    static const std::size_t F = 2;
    static const std::size_t G = 3;
    static const std::size_t P = 4;
    static const std::size_t NUM_VARIABLES = 5;

    typedef hpx::serialization::serialize_buffer<Real> buffer_type;
    typedef std::vector<hpx::shared_future<void> > future_vector;
    typedef std::vector<future_vector> future_grid;


    partition_server() {}
    ~partition_server() {}

    partition_server(partition_server const& other)
    : c(other.c)
    {
        std::cout << "constr" << std::endl;
    }

    partition_server(partition_server && other)
    : c(other.c),
      is_left_(other.is_left_),
      is_right_(other.is_right_),
      is_bottom_(other.is_bottom_),
      is_top_(other.is_top_),
      current(other.current),
      last(other.last),
      rhs_data_(std::move(other.rhs_data_)),
      cell_type_data_(other.cell_type_data_),
      fluid_cells_(other.fluid_cells_),
      boundary_cells_(other.boundary_cells_),
      obstacle_cells_(other.obstacle_cells_),
      cells_x_(other.cells_x_),
      cells_y_(other.cells_y_),
      idx_(other.idx_),
      idy_(other.idy_),
      step_(other.step_),
      outcount_(other.outcount_),
      t_(other.t_),
      next_out_(other.next_out_)
    {
        for (std::size_t var = 0; var < NUM_VARIABLES; ++var)
            data_[var] = std::move(other.data_[var]);
    }

    partition_server(io::config const& cfg);

    void init();
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, init, init_action);

    /// return the sliced data appropriate for given direction
    std::pair<Real, Real> do_timestep(Real dt);
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, do_timestep, do_timestep_action);

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

    send_buffer<buffer_type, BOTTOM_RIGHT, set_top_left_boundary_action> send_buffer_bottom_right_;
    recv_buffer<buffer_type, BOTTOM_RIGHT> recv_buffer_bottom_right_[NUM_VARIABLES];

    send_buffer<buffer_type, TOP_LEFT, set_bottom_right_boundary_action> send_buffer_top_left_;
    recv_buffer<buffer_type, TOP_LEFT> recv_buffer_top_left_[NUM_VARIABLES];

protected:
    template<direction dir>
    void send_boundary(std::size_t step, std::size_t var, future_vector& send_future);

    template<std::size_t var>
    void send_right_and_top_boundaries(std::size_t step, future_grid& send_futures);

    template<std::size_t var>
    void send_cross_boundaries(std::size_t step, future_grid& send_futures);

    template<std::size_t var>
    void send_all_boundaries(std::size_t step, future_grid& send_futures);

    template<direction dir>
    void receive_boundary(std::size_t step, std::size_t var, future_vector& recv_futures);

    template<std::size_t var>
    void receive_left_and_bottom_boundaries(std::size_t step, future_grid& recv_futures);

    template<std::size_t var>
    void receive_cross_boundaries(std::size_t step, future_grid& recv_futures);

    template<std::size_t var>
    void receive_all_boundaries(std::size_t step, future_grid& recv_futures);

    void wait_all_boundaries(future_grid& recv_futures);

    template<direction dir>
    hpx::shared_future<void> get_dependency(std::size_t idx_block,
        std::size_t idy_block, future_grid const& recv_futures,
        partition_data<hpx::shared_future<void> > const& calc_futures);

    void init_future_grid(future_grid& grid)
    {
        grid[LEFT].resize(c.num_y_blocks);
        grid[RIGHT].resize(c.num_y_blocks);
        grid[BOTTOM].resize(c.num_x_blocks);
        grid[TOP].resize(c.num_x_blocks);
        grid[BOTTOM_RIGHT].resize(1);
        grid[TOP_LEFT].resize(1);
    }

    void clear_and_reserve(future_grid& grid)
    {
        for (std::size_t dir = 0; dir != NUM_DIRECTIONS; ++dir)
            grid[dir].clear();

        grid[LEFT].reserve(c.num_y_blocks);
        grid[RIGHT].reserve(c.num_y_blocks);
        grid[BOTTOM].reserve(c.num_x_blocks);
        grid[TOP].reserve(c.num_x_blocks);
        grid[BOTTOM_RIGHT].reserve(1);
        grid[TOP_LEFT].reserve(1);
    }

    void clear(future_grid& grid)
    {
        for (std::size_t dir = 0; dir != NUM_DIRECTIONS; ++dir)
            grid[dir].clear();
    }

private:

    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned version)
    {
        ar & c & is_left_ & is_right_ & is_bottom_ & is_top_ & current & last & data_ & rhs_data_ & cell_type_data_
           & fluid_cells_ & boundary_cells_ & obstacle_cells_ & cells_x_ & cells_y_ & idx_ & idy_
           & step_ & outcount_ & t_ & next_out_;
    }

    partition_data<Real> data_[NUM_VARIABLES];
    partition_data<Real> rhs_data_;
    partition_data<std::bitset<6> > cell_type_data_;

    partition_data<std::vector<std::pair<std::size_t, std::size_t> > > fluid_cells_;
    partition_data<std::vector<std::pair<std::size_t, std::size_t> > > boundary_cells_;
    partition_data<std::vector<std::pair<std::size_t, std::size_t> > > obstacle_cells_;

    partition_data<hpx::shared_future<void> > set_velocity_futures;
    future_grid send_futures_U;
    future_grid send_futures_V;
    future_grid recv_futures_U;
    future_grid recv_futures_V;

    partition_data<hpx::shared_future<void> > compute_fg_futures;
    future_grid send_futures_F;
    future_grid send_futures_G;
    future_grid recv_futures_F;
    future_grid recv_futures_G;

    partition_data<hpx::shared_future<void> > compute_rhs_futures;

    partition_data<hpx::shared_future<Real> > compute_res_futures;

    std::size_t current;
    std::size_t last;

    partition_data<hpx::shared_future<void> > set_p_futures;
    future_grid send_futures_P;
    std::array<future_grid, 2> recv_futures_P;
    std::array<partition_data<hpx::shared_future<void> >, 2> sor_cycle_futures;

    std::vector<hpx::future<std::pair<Real, Real> > > local_max_uvs;

    std::vector<hpx::id_type> ids_;

    io::config c;

    std::size_t cells_x_, cells_y_;
    std::size_t idx_, idy_;
    std::size_t step_;
    std::size_t outcount_;

    Real t_, next_out_;

    util::cancellation_token token;

    bool is_left_, is_right_, is_bottom_, is_top_;
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
