#ifndef COMPUTATION_CUSTOM_GRAIN_SIZE_HPP
#define COMPUTATION_CUSTOM_GRAIN_SIZE_HPP

#include "strategy.hpp"

namespace computation {

class custom_grain_size : public strategy {

    public:
        using strategy::strategy;
        virtual void set_boundary(vector_grid_type& uv_grid);
        virtual void set_pressure_on_boundary(scalar_grid_type& p_grid);
        virtual void compute_fg(vector_grid_type& fg_grid, vector_grid_type const& uv_grid, RealType dt);
        virtual void compute_rhs(scalar_grid_type& rhs_grid, vector_grid_type const& fg_grid, RealType dt);

        virtual void sor_cycle(scalar_grid_type& p_grid, scalar_grid_type const& rhs_grid);
        virtual hpx::future<RealType> compute_residual(scalar_grid_type const& p_grid, scalar_grid_type const& rhs_grid);

        virtual hpx::future<std::pair<RealType, RealType> > update_velocities(vector_grid_type& uv_grid,
                                                        vector_grid_type const& fg_grid, scalar_grid_type const& p_grid, RealType dt);
};

}//namespace

#endif
