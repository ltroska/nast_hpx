#ifndef NAST_HPX_GRID_SERVER_PARTITION_SERVER_HPP
#define NAST_HPX_GRID_SERVER_PARTITION_SERVER_HPP

#include <hpx/include/components.hpp>
#include <hpx/runtime/serialization/serialize.hpp>

#include "grid/partition_data.hpp"
#include "grid/send_buffer.hpp"
#include "grid/recv_buffer.hpp"

#include "io/config.hpp"

namespace nast_hpx { namespace grid { namespace server {

char const* partition_basename = "/nast_hpx/partition/";
char const* residual_basename = "/nast/hpx/partition/residual";

    
/// component encapsulates partition_data, making it remotely available
struct HPX_COMPONENT_EXPORT partition_server
    : public hpx::components::component_base<partition_server>
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

    partition_server(io::config&& cfg, uint idx, uint idy);
        
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
    void receive_left_and_bottom_boundaries(std::size_t step, future_vector& recv_futures);
     
    template<std::size_t var>
    void receive_cross_boundaries(std::size_t step, future_vector& recv_futures);
    
    template<std::size_t var>
    void receive_all_boundaries(std::size_t step, future_vector& recv_futures);
    
    void wait_all_boundaries(future_vector& recv_futures);
    
    template<direction dir>
    hpx::shared_future<void> get_dependency(std::size_t idx_block,
        std::size_t idy_block, future_vector const& recv_futures,
        partition_data<hpx::shared_future<void> > const& calc_futures);

private:
    
    partition_data<Real> data_[NUM_VARIABLES];
    partition_data<Real> rhs_data_;
    partition_data<std::bitset<6> > cell_types_;
    std::vector<hpx::id_type> ids_;
    std::size_t rank_, num_partitions_, num_partitions_x_, num_partitions_y_;
    std::size_t cells_x_, cells_y_;
    std::size_t num_x_blocks_, num_y_blocks_, cells_per_x_block_, cells_per_y_block_;
    std::size_t num_fluid_cells_;
    std::size_t idx_, idy_;
    std::size_t step_;
    std::size_t i_max_, j_max_;
    std::size_t iter_max_;
    std::size_t outcount_;
        
    Real dx_, dy_, dx_sq_, dy_sq_, over_dx_sq_, over_dy_sq_, over_dx_, over_dy_,
        part1_, part2_, alpha_, beta_, re_, gx_, gy_, eps_sq_, t_, delta_vec_,
        next_out_;
    
    boundary_data u_bnd_, v_bnd_, uv_bnd_type_;
    
    bool is_left_, is_right_, is_bottom_, is_top_, vtk_;
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
