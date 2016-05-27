#ifndef NAST_HPX_STEPPER_STEPPER_HPP
#define NAST_HPX_STEPPER_STEPPER_HPP

#include "server/stepper_server.hpp"

namespace nast_hpx { namespace stepper {
    
/// A client that represents the stepper_server component.

struct stepper
    : hpx::components::client_base<stepper, server::stepper_server>
{
    typedef hpx::components::client_base<stepper, server::stepper_server>
        base_type;

    // Contruct a new stepper client.
    stepper()
      : base_type(hpx::new_<server::stepper_server>
          (hpx::find_here(), hpx::get_num_localities_sync()))
    {
        hpx::register_with_basename(server::stepper_basename, get_id(),
                                        hpx::get_locality_id());
    }

    // Construct new instances/wrap existing steppers from other localities.
    stepper(hpx::id_type loc)
      : base_type(hpx::new_<server::stepper_server>
          (loc, hpx::get_num_localities_sync()))
    {
        hpx::register_with_basename(
            server::stepper_basename, get_id(), hpx::get_locality_id());
    }

    // Method to wrap existing steppers, since a stepper is component is a 
    // future.
    stepper(hpx::future<hpx::id_type> && id)
      : base_type(std::move(id))
    {}

    // Method forwards the config to the wrapped stepper_server.
    hpx::future<void> setup(io::config&& cfg)
    {
        server::stepper_server::setup_action act;
        return hpx::async(act, get_id(), std::move(cfg));
    }
};


}//namespace stepper
}

#endif
