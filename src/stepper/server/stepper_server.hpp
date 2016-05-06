#ifndef STEPPER_SERVER_STEPPER_HPP
#define STEPPER_SERVER_STEPPER_HPP

#include <hpx/include/components.hpp>
#include <hpx/lcos/local/receive_buffer.hpp>

#include "io/config.hpp"

#ifdef CUSTOM_GRAIN_SIZE
#include "computation/custom_grain_size.hpp"
#else
#include "computation/with_for_each.hpp"
#endif

#include "grid/types.hpp"

#include <hpx/error.hpp>

namespace stepper { namespace server {

char const* stepper_basename = "/nast_hpx/stepper/";
char const* residual_basename = "/nast_hpx/gather/residual";
char const* velocity_basename = "/nast_hpx/gather/velocity";
char const* dt_basename = "/nast_hpx/gather/dt";

/// Component responsible for the timestepping and communication of data.
struct HPX_COMPONENT_EXPORT stepper_server
    : hpx::components::component_base<stepper_server>
{
    public:
        stepper_server() {}
        stepper_server(uint num_localities);

        /// Method that sets up the stepper with a given config
        void setup(io::config&& cfg);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, setup, setup_action);

        /// Method that starts the time stepping loop
        void do_work();

        /// Method that does a single timestep
        hpx::future<std::pair<RealType, RealType> > do_timestep(uint step,
                                                                RealType dt);

        /// Action that places the bool into a receive buffer.
        /// In iteration n the SOR loop takes the n-th value from the buffer
        /// and decides to start another iteration / stops, depending on the
        /// value.
        void set_keep_running(uint iter, bool kr);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, set_keep_running,
                                        set_keep_running_action);        
                                        
        void set_dt(uint step, RealType dt);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, set_dt, set_dt_action);

        /// Action that places a component representing the pressure block of a
        /// neighboring locality into a buffer. The data was sent from that
        /// neighboring locality into given direction.
        void receive_p_action_(uint t, scalar_partition p, direction to_dir);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, receive_p_action_,
                                        receive_p_action);

        /// Action that places a component representing the FG value block of a
        /// neighboring locality into a buffer. The data was sent from that
        /// neighboring locality into given direction.
        void receive_fg_action_(uint t, vector_partition fg, direction to_dir);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, receive_fg_action_,
                                        receive_fg_action);

        /// Helper action that places a component representing the velocity
        /// block of a neighboring locality into a buffer. The data was sent 
        /// from that neighboring locality into given direction.
        void receive_uv_action_(uint t, vector_partition uv, direction to_dir);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, receive_uv_action_,
                                        receive_uv_action);

    protected:
        
        /// Method that prints the given grid asynchronously on std::cout.
        template<typename T>
        void print_grid(std::vector<grid::partition<T> > const& grid,
                std::string const& message = "") const;
        
        /// Method that writes all grids to a vtk file.
        void write_vtk(uint step);

        /// Method that converts 2D to 1D indices.
        uint get_index(uint k, uint l) const;

        /// Helper action that sends a pressure block to the locality in the 
        /// given direction.
        void send_p_to_neighbor(uint t, scalar_partition p, direction dir);
        /// Helper method that retrieves pressure blocks from the buffer
        /// corresponding to the locality in the given direction.
        scalar_partition receive_p_from_neighbor(uint t, direction dir);
        /// Method that manages the communication of all pressure data
        /// (sending and receiving) for the given step.
        void communicate_p_grid(uint iter);
        
        /// Helper action that sends a FG block to the locality in the given
        /// direction.
        void send_fg_to_neighbor(uint t, vector_partition fg, direction dir);
        /// Helper method that retrieves FG blocks from the buffer
        /// corresponding to the locality in the given direction.
        vector_partition receive_fg_from_neighbor(uint t, direction dir);
        /// Method that manages the communication of all FG data
        /// (sending and receiving) for the given step.
        void communicate_fg_grid(uint step);

        /// Helper action that sends a velocity block to the locality in the 
        /// given direction.
        void send_uv_to_neighbor(uint t, vector_partition uv, direction dir);
        /// Helper method that retrieves velocity blocks from the buffer
        /// corresponding to the locality in the given direction.
        vector_partition receive_uv_from_neighbor(uint t, direction dir);
        /// Method that manages the communication of all velocity data
        /// (sending and receiving) for the given step.
        void communicate_uv_grid(uint step);

    private:
        void initialize_parameters();
        void initialize_grids();
        void initialize_communication();

        io::config c;
        
        struct {
            uint i_max, j_max, num_partitions_x, num_partitions_y,
                num_cells_per_partition_x, num_cells_per_partition_y;
            
            RealType re, pr, omega, alpha, beta, dx, dy;
        } params;
        
        // we decide on compile time which computation strategy to use
#ifdef CUSTOM_GRAIN_SIZE
        typedef computation::custom_grain_size strategy;
#else
        typedef computation::with_for_each strategy;
#endif

        uint num_localities, num_localities_x, num_localities_y;
        std::vector<hpx::naming::id_type> localities;

        // grid representing the beginning global indices of each block
        index_grid_type index_grid;
        
        // velocity and FG grid
        vector_grid_type uv_grid, fg_grid;
        
        // temporary grids used in the computation
        vector_grid_type uv_temp_grid, fg_temp_grid;

        // scalar grids
        scalar_grid_type p_grid, rhs_grid, temperature_grid, stream_grid,
            vorticity_grid, heat_grid;
        
        // corresponding needed temporary grids
        scalar_grid_type p_temp_grid, temperature_temp_grid;

        // grid describing the type of each cell (boundary, obstacle, fluid)
        flag_grid_type flag_grid;
        
        std::vector<std::vector<std::pair<uint, uint> > > obstacle;
        std::vector<std::vector<std::vector<std::pair<uint, uint> > > > boundary;
        std::vector<std::vector<std::pair<uint, uint> > > fluid;

        RealType t, next_write;
        uint out_iter;

        // arrays for neighbor localities
        bool has_neighbor[NUM_DIRECTIONS];
        hpx::shared_future<hpx::id_type> neighbor_steppers_[NUM_DIRECTIONS];
        
        // buffer for the data of the neighboring localities
        hpx::lcos::local::receive_buffer<scalar_partition>
            p_recv_buffs_[NUM_DIRECTIONS];
        
        hpx::lcos::local::receive_buffer<vector_partition>
            fg_recv_buffs_[NUM_DIRECTIONS];
        
        hpx::lcos::local::receive_buffer<vector_partition>
            uv_recv_buffs_[NUM_DIRECTIONS];
        
        // buffer for the SOR loop
        hpx::lcos::local::receive_buffer<bool> keep_running;
        
        hpx::lcos::local::receive_buffer<RealType> dt_buffer;

        // dummies so we have a parameter for dataflow if this locality does
        // not have full amount of neighbors
        scalar_partition scalar_dummy;
        vector_partition vector_dummy;
};

}//namespace server
}//namespace stepper

HPX_REGISTER_ACTION_DECLARATION(stepper::server::stepper_server::setup_action,
                                    stepper_server_setup_action);

#endif
