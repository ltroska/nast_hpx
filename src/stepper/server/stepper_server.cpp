#include <cmath>
#include <chrono>

#include <hpx/include/iostreams.hpp>

#include "stepper_server.hpp"

typedef nast_hpx::stepper::server::stepper_server stepper_component;
typedef hpx::components::component<stepper_component> stepper_server_type;

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_COMPONENT(stepper_server_type, stepper_component);
HPX_REGISTER_ACTION(nast_hpx::stepper::server::stepper_server::setup_action,
    stepper_server_setup_action);

namespace nast_hpx { namespace stepper { namespace server {

stepper_server::stepper_server(uint nl)
: num_localities(nl)
{}

void stepper_server::setup(io::config&& cfg)
{
    // special case for two localities, we want a square configuration
    // of localities
    if (num_localities == 2)
    {
        num_localities_x = 2;
        num_localities_y = 1;
    }
    else
    {
        num_localities_x = static_cast<uint> (sqrt(num_localities));
        num_localities_y = num_localities_x;
    }


    std::cout << "hi" << std::endl;
    part = grid::partition(hpx::find_here(), 16, 16);
    part.do_timestep();

   
}


}//namespace server
}//namespace stepper
}