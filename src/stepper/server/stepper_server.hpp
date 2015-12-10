#ifndef STEPPER_SERVER_STEPPER_HPP
#define STEPPER_SERVER_STEPPER_HPP

#include <hpx/include/components.hpp>
#include <hpx/hpx.hpp>

#include "grid/partition.hpp"

namespace stepper { namespace server {

char const* stepper_basename = "/cfd_hpx/stepper/";

enum direction
{
    left, right, top, bottom, top_left, top_right, bottom_left, bottom_right
};

enum stepper_type
{
    left_boundary, right_boundary, top_boundary, bottom_boundary, top_left_boundary, top_right_boundary, bottom_left_boundary, bottom_right_boundary, interior
};

struct HPX_COMPONENT_EXPORT stepper_server
    : hpx::components::component_base<stepper_server>
{
    public:
        typedef std::vector<std::vector<grid::partition> > space;

        stepper_server() {}
        stepper_server(uint num_localities);

        uint do_work(uint num_local_partitions_x, uint num_local_partitions_y, uint num_cells_x, uint num_cells_y, RealType delta_x, RealType delta_y);

        HPX_DEFINE_COMPONENT_ACTION(stepper_server, do_work, do_work_action);

        uint setup(uint i_max, uint j_max, RealType x_length, RealType y_length,  uint num_partitions_x, uint num_partitions_y);

        HPX_DEFINE_COMPONENT_ACTION(stepper_server, setup, setup_action);

        void receive_from_neighbor(uint t, grid::partition p, direction dir)
        {
            switch (dir)
            {
                case left:
                    U[num_local_partitions_x_ - 1][t] = std::move(p);
                    break;
                case right:
                    U[0][t] = std::move(p);
                    break;
                case top:
                    U[t][0] = std::move(p);
                    break;
                case bottom:
                    U[t][num_local_partitions_y_ - 1] = std::move(p);
                    break;
                case top_left:
                    U[num_local_partitions_x_ - 1][0] = std::move(p);
                    break;
                case top_right:
                    U[0][0] = std::move(p);
                    break;
                case bottom_left:
                    U[num_local_partitions_x_ - 1][num_local_partitions_y_ - 1] = std::move(p);
                    break;
                case bottom_right:
                    U[0][num_local_partitions_y_ - 1] = std::move(p);
                    break;
                default:
                    return;
            }
        }

        HPX_DEFINE_COMPONENT_ACTION(stepper_server, receive_from_neighbor, receive_action);

    protected:

        void set_boundary_values_u_v();

        void write_vtk_files();

        void send_to_neighbor(uint t, grid::partition p, direction dir)
        {
            hpx::id_type neighbor;

            switch (dir)
            {
                case left:
                    neighbor = left_.get();
                    break;
                case right:
                    neighbor = right_.get();
                    break;
                case top:
                    neighbor = top_.get();
                    break;
                case bottom:
                    neighbor = bottom_.get();
                    break;
                case top_left:
                    neighbor = top_left_.get();
                    break;
                case top_right:
                    neighbor = top_right_.get();
                    break;
                case bottom_left:
                    neighbor = bottom_left_.get();
                    break;
                case bottom_right:
                    neighbor = bottom_right_.get();
                    break;
                default:
                    return;
            }

            hpx::apply(receive_action(), neighbor, t, p, dir);
        }

    private:
        space U;

        RealType dx_;
        RealType dy_;

        RealType x_length_;
        RealType y_length_;

        uint num_local_partitions_x_;
        uint num_local_partitions_y_;
        uint num_cells_x_;
        uint num_cells_y_;
        uint res_x_;
        uint res_y_;

        uint num_localities_;
        uint locality_id_;
        stepper_type type_;

        hpx::shared_future<hpx::id_type> left_, right_, top_, bottom_, top_left_, top_right_, bottom_left_, bottom_right_;
};

}//namespace server
}//namespace stepper

HPX_REGISTER_ACTION_DECLARATION(stepper::server::stepper_server::do_work_action, stepper_server_do_work_action);
HPX_REGISTER_ACTION_DECLARATION(stepper::server::stepper_server::setup_action, stepper_server_setup_action);
#endif
