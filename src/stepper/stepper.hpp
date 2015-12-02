#ifndef STEPPER_STEPPER_HPP
#define STEPPER_STEPPER_HPP

#include "server/stepper_server.hpp"


namespace stepper {

struct stepper
    : hpx::components::client_base<stepper, server::stepper_server>
{
    typedef hpx::components::client_base<stepper, server::stepper_server> base_type;

    // construct new instances/wrap existing steppers from other localities
    stepper()
      : base_type(hpx::new_<server::stepper_server>
          (hpx::find_here()))
    {
       // hpx::register_with_basename(
        //    stepper_basename, get_id(), hpx::get_locality_id());
    }

    stepper(hpx::future<hpx::id_type> && id)
      : base_type(std::move(id))
    {}

    hpx::future<server::stepper_server::space> do_work()
    {
        server::stepper_server::do_work_action act;
        return hpx::async(act, get_id());
    }
};


}//namespace stepper

#endif
