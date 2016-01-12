#ifndef STEPPER_STEPPER_HPP
#define STEPPER_STEPPER_HPP

#include "server/stepper_server.hpp"

namespace stepper {

struct stepper
    : hpx::components::client_base<stepper, server::stepper_server>
{
    typedef hpx::components::client_base<stepper, server::stepper_server> base_type;

    stepper()
      : base_type(hpx::new_<server::stepper_server>
          (hpx::find_here(), hpx::get_num_localities_sync()))
    {
        hpx::register_with_basename(
            server::stepper_basename, get_id(), hpx::get_locality_id());
    }

    // construct new instances/wrap existing steppers from other localities
    stepper(hpx::id_type loc)
      : base_type(hpx::new_<server::stepper_server>
          (loc, hpx::get_num_localities_sync()))
    {
        hpx::register_with_basename(
            server::stepper_basename, get_id(), hpx::get_locality_id());
    }

    stepper(hpx::future<hpx::id_type> && id)
      : base_type(std::move(id))
    {}

    hpx::future<uint> setup(uint i_max, uint j_max, RealType x_length, RealType y_length, uint num_partitions_x, uint num_partitions_y)
    {
        server::stepper_server::setup_action act;
        return hpx::async(act, get_id(), i_max, j_max, x_length, y_length, num_partitions_x, num_partitions_y);
    }

    hpx::future<uint> setup(cfd_config config)
    {
        server::stepper_server::setup_with_config_action act;
        return hpx::async(act, get_id(), config);
    }

    hpx::future<uint> update_delta_t(RealType dt)
    {
        server::stepper_server::update_delta_t_action act;
        return hpx::async(act, get_id(), dt);
    }

    hpx::future<uint> set_velocity_on_boundary()
    {
        server::stepper_server::set_velocity_action act;
        return hpx::async(act, get_id());
    }

    hpx::future<uint> set_pressure_on_boundary()
    {
        server::stepper_server::set_pressure_action act;
        return hpx::async(act, get_id());
    }

    hpx::future<uint> compute_fg()
    {
        server::stepper_server::compute_fg_action act;
        return hpx::async(act, get_id());
    }

    hpx::future<uint> set_rhs()
    {
        server::stepper_server::set_rhs_action act;
        return hpx::async(act, get_id());
    }

    hpx::future<RealType> get_residual()
    {
        server::stepper_server::get_residual_action act;
        return hpx::async(act, get_id());
    }

    hpx::future<grid::vector_cell> update_velocities()
    {
        server::stepper_server::update_velocities_action act;
        return hpx::async(act, get_id());
    }

    hpx::future<uint> sor_cycle()
    {
        server::stepper_server::sor_cycle_action act;
        return hpx::async(act, get_id());
    }

    hpx::future<uint> do_work(uint num_local_partitions_x, uint num_local_partitions_y, uint num_cells_x, uint num_cells_y, RealType dx, RealType dy)
    {
        server::stepper_server::do_work_action act;
        return hpx::async(act, get_id(), num_local_partitions_x, num_local_partitions_y, num_cells_x, num_cells_y, dx, dy);
    }

};


}//namespace stepper

#endif
