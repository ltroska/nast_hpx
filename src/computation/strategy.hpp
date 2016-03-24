#ifndef COMPUTATION_STRATEGY_HPP
#define COMPUTATION_STRATEGY_HPP

#include <hpx/hpx.hpp>

#include "grid/types.hpp"

namespace computation {

class strategy {
    public:
        virtual void set_velocity_on_boundary(vector_data& uv_data, uint global_i, uint global_j,
                                                uint i_max, uint j_max) {}

};

}//namespace

#endif
