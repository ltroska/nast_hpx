#ifndef NAST_HPX_STEPPER_SERVER_STEPPER_HPP
#define NAST_HPX_STEPPER_SERVER_STEPPER_HPP

#include <hpx/include/components.hpp>
#include <hpx/error.hpp>

#include "io/config.hpp"
#include "grid/partition.hpp"

namespace nast_hpx { namespace stepper { namespace server {

char const* stepper_basename = "/nast_hpx/stepper/";
char const* velocity_basename = "/nast_hpx/gather/velocity";
char const* barrier_basename = "/nast_hpx/barrier";

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

        void set_dt(uint step, Real dt);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, set_dt, set_dt_action);

        void do_timestep(Real dt);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, do_timestep, do_timestep_action);

    private:
        uint num_localities, num_localities_x, num_localities_y;
        hpx::lcos::local::receive_buffer<Real> dt_buffer;
};

}//namespace server
}//namespace stepper
}

HPX_REGISTER_ACTION_DECLARATION(nast_hpx::stepper::server::stepper_server::setup_action,
                                    stepper_server_setup_action);

#endif
