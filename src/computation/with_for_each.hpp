#ifndef WITH_FOR_EACH_HPP
#define WITH_FOR_EACH_HPP

#include "strategy.hpp"

namespace computation {

class with_for_each : public strategy {

    public:
        virtual void set_velocity_on_boundary(vector_data& uv_data, uint global_i, uint global_j,
                                                uint i_max, uint j_max);
};

}//namespace

#endif
