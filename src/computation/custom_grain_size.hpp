#ifndef COMPUTATION_CUSTOM_GRAIN_SIZE_HPP
#define COMPUTATION_CUSTOM_GRAIN_SIZE_HPP

#include "strategy.hpp"

namespace computation {

class custom_grain_size : public strategy {

    public:
        using strategy::strategy;
        virtual void set_velocity_on_boundary(vector_grid_type& uv_grid);
        virtual void compute_fg(vector_grid_type& fg_grid, vector_grid_type const& uv_grid, RealType dt);

};

}//namespace

#endif
