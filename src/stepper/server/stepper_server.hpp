#ifndef STEPPER_SERVER_STEPPER_HPP
#define STEPPER_SERVER_STEPPER_HPP

#include <hpx/include/components.hpp>
#include <hpx/hpx.hpp>

#include "grid/partition.hpp"

namespace stepper { namespace server {

char const* stepper_basename = "/cfd_hpx/stepper/";

struct HPX_COMPONENT_EXPORT stepper_server
    : hpx::components::component_base<stepper_server>
{
    public:
        typedef std::vector<std::vector<grid::partition> > space;

        stepper_server() {}
        stepper_server(uint num_localities);

        uint do_work(uint num_local_partitions_x, uint num_local_partitions_y, uint num_cells_x, uint num_cells_y, RealType delta_x, RealType delta_y);

        HPX_DEFINE_COMPONENT_ACTION(stepper_server, do_work, do_work_action);

        void from_right(uint t, grid::partition p)
        {
            right_receive_buffer_.store_received(t, std::move(p));
        }

        HPX_DEFINE_COMPONENT_ACTION(stepper_server, from_right, from_right_action);

    protected:

        void set_boundary_values_u_v();

        void write_vtk_files();

        grid::partition receive_right(uint t)
        {
            return right_receive_buffer_.receive(t).get();
        }

        void send_left(uint t, grid::partition p)
        {
            hpx::apply(from_right_action(), left_.get(), t, p);
            hpx::cout << "done senidng" << hpx::endl << hpx::flush;
        }

    private:
        space U;

        RealType dx;
        RealType dy;

        uint num_local_partitions_x_;
        uint num_local_partitions_y_;
        uint num_cells_x_;
        uint num_cells_y_;

        uint num_localities_;

        hpx::shared_future<hpx::id_type> left_;
        hpx::lcos::local::receive_buffer<grid::partition> right_receive_buffer_;
};

}//namespace server
}//namespace stepper

HPX_REGISTER_ACTION_DECLARATION(stepper::server::stepper_server::do_work_action, stepper_server_do_work_action);
#endif
