#ifndef WITH_FOR_EACH_HPP
#define WITH_FOR_EACH_HPP

#include "strategy.hpp"

namespace computation {

class with_for_each : public strategy {

    public:
        virtual void set_velocity_on_boundary(vector_data& uv, uint global_i, uint global_j,
                                                uint i_max, uint j_max);

        virtual void compute_fg(vector_data& fg, vector_data const& uv_center,
                            vector_data const& uv_left, vector_data const& uv_right,
                            vector_data const& uv_bottom, vector_data const& uv_top,
                            vector_data const& uv_bottomright, vector_data const& uv_topleft,
                            uint global_i, uint global_j, uint i_max, uint j_max, RealType re,
                            RealType dx, RealType dy, RealType dt, RealType alpha);

        virtual void compute_rhs(scalar_data& rhs, vector_data const& fg_center, vector_data const& fg_left,
                                    vector_data const& fg_bottom, uint global_i, uint global_j, uint i_max, uint j_max,
                                    RealType dx, RealType dy, RealType dt);
};

}//namespace

#endif
