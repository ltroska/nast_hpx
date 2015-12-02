#ifndef STEPPER_SERVER_STEPPER_HPP
#define STEPPER_SERVER_STEPPER_HPP

#include <hpx/include/components.hpp>

#include "grid/partition.hpp"

namespace stepper { namespace server {

struct HPX_COMPONENT_EXPORT stepper_server
    : hpx::components::component_base<stepper_server>
{
    public:
        typedef std::vector<std::vector<grid::partition> > space;

        stepper_server();

        stepper_server(uint num_local_partitions_x, uint num_local_partitions_y, uint num_cells_x, uint num_cells_y, RealType delta_x, RealType delta_y);

        uint do_work();

        HPX_DEFINE_COMPONENT_ACTION(stepper_server, do_work, do_work_action);

    protected:

        void set_boundary_values_u_v();

        void write_vtk_files();

    private:
        space U;

        RealType dx;
        RealType dy;

        uint num_local_partitions_x;
        uint num_local_partitions_y;
        uint num_cells_x;
        uint num_cells_y;
};

}//namespace server
}//namespace stepper

HPX_REGISTER_ACTION_DECLARATION(stepper::server::stepper_server::do_work_action, stepper_server_do_work_action);
#endif
