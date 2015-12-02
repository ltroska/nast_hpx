#ifndef STEPPER_SERVER_STEPPER_HPP
#define STEPPER_SERVER_STEPPER_HPP

#include <hpx/include/components.hpp>

#include "grid/partition.hpp"
#include "internal/cfd_config.hpp"

namespace stepper { namespace server {

struct HPX_COMPONENT_EXPORT stepper_server
    : hpx::components::component_base<stepper_server>
{
    public:
        typedef std::vector<std::vector<grid::partition> > space;

        stepper_server();

        space do_work();

        HPX_DEFINE_COMPONENT_ACTION(stepper_server, do_work, do_work_action);

    protected:

        void write_vtk_file() {}

    private:
        space U;

        RealType dx;
        RealType dy;

        uint numPartitions;
        uint numLocalities;


};

}//namespace server
}//namespace stepper

HPX_REGISTER_ACTION_DECLARATION(stepper::server::stepper_server::do_work_action, stepper_server_do_work_action);
#endif
