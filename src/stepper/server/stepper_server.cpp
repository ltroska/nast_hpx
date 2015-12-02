#include "stepper_server.hpp"

typedef stepper::server::stepper_server stepper_component;
typedef hpx::components::component<stepper_component> stepper_server_type;

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_COMPONENT(stepper_server_type, stepper_component);

HPX_REGISTER_ACTION(stepper::server::stepper_server::do_work_action, stepper_server_do_work_action);


namespace stepper { namespace server {
//calculate locality to distribute partitions evenly
inline uint locidx(uint i, uint np, uint nl)
{
    return i / (np/nl);
}

stepper_server::stepper_server()
{
   // hpx::cout << hpx::find_here() << hpx::flush;
}

stepper_server::space stepper_server::do_work()
{
   // hpx::cout << hpx::find_here() << hpx::flush()
}

}//namespace server
}//namespace stepper
