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


    auto rank = hpx::get_locality_id();

    std::cout << cfg << std::endl;
    std::size_t idx, idy;
    
    
    idx = (rank % cfg.num_localities_x);
    idy = (rank / cfg.num_localities_x);
    part = grid::partition(hpx::find_here(), cfg, idx, idy);
    part.init_sync();
    part.do_timestep(0.001);

   
}


}//namespace server
}//namespace stepper
}