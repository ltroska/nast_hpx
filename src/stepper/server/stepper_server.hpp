#ifndef STEPPER_SERVER_STEPPER_HPP
#define STEPPER_SERVER_STEPPER_HPP

#include <hpx/include/components.hpp>
#include <hpx/lcos/local/receive_buffer.hpp>

#include "io/config.hpp"
#include "computation/parameters.hpp"
#include "computation/strategy.hpp"
#include "grid/types.hpp"
#include "util/helpers.hpp"

#include <hpx/error.hpp>

namespace stepper { namespace server {

char const* stepper_basename = "/nast_hpx/stepper/";
char const* gather_basename = "/nast_hpx/gather/";

struct HPX_COMPONENT_EXPORT stepper_server
    : hpx::components::component_base<stepper_server>
{
    public:
        stepper_server() {}
        stepper_server(uint num_localities);

        void setup(io::config cfg);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, setup, setup_action);

        void do_work();

        std::pair<RealType, RealType> do_timestep(uint step, RealType dt);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, do_timestep, do_timestep_action);

        void do_sor_cycle();
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, do_sor_cycle, do_sor_cycle_action);

        void set_keep_running(uint iter, bool kr);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, set_keep_running, set_keep_running_action);

        void receive_p_action_(uint t, scalar_partition p, direction to_dir);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, receive_p_action_, receive_p_action);

        void receive_fg_action_(uint t, vector_partition fg, direction to_dir);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, receive_fg_action_, receive_fg_action);

        void receive_uv_action_(uint t, vector_partition uv, direction to_dir);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, receive_uv_action_, receive_uv_action);

    protected:
        template<typename T>
        void print_grid(std::vector<grid::partition<T> > const& grid, const std::string message = "") const;
        void write_vtk(uint step);

        uint get_index(uint k, uint l) const;

        void sor();

        void send_p_to_neighbor(uint t, scalar_partition p, direction dir);
        scalar_partition receive_p_from_neighbor(uint t, direction dir);
        void communicate_p_grid(uint iter);

        void send_fg_to_neighbor(uint t, vector_partition fg, direction dir);
        vector_partition receive_fg_from_neighbor(uint t, direction dir);
        void communicate_fg_grid(uint step);

        void send_uv_to_neighbor(uint t, vector_partition uv, direction dir);
        vector_partition receive_uv_from_neighbor(uint t, direction dir);
        void communicate_uv_grid(uint step);

    private:
        void initialize_parameters();
        void initialize_grids();
        void initialize_communication();

        RealType compute_new_dt(std::pair<RealType, RealType>) const;

        io::config c;
        computation::strategy* strategy;
        computation::parameters params;

        uint num_localities, num_localities_x, num_localities_y;
        std::vector<hpx::naming::id_type> localities;

        index_grid_type index_grid;
        vector_grid_type uv_grid, fg_grid;
        scalar_grid_type p_grid, rhs_grid, temperature_grid, stream_grid, vorticity_grid, heat_grid;
        flag_grid_type flag_grid;

        RealType t, next_write;
        uint out_iter;

        bool has_neighbor[NUM_DIRECTIONS];
        hpx::shared_future<hpx::id_type> neighbor_steppers_[NUM_DIRECTIONS];
        hpx::lcos::local::receive_buffer<scalar_partition> p_recv_buffs_[NUM_DIRECTIONS];
        hpx::lcos::local::receive_buffer<vector_partition> fg_recv_buffs_[NUM_DIRECTIONS];
        hpx::lcos::local::receive_buffer<vector_partition> uv_recv_buffs_[NUM_DIRECTIONS];
        hpx::lcos::local::receive_buffer<bool> keep_running;

        scalar_partition scalar_dummy;
        vector_partition vector_dummy;
};

}//namespace server
}//namespace stepper

HPX_REGISTER_ACTION_DECLARATION(stepper::server::stepper_server::setup_action, stepper_server_setup_action);

#endif
