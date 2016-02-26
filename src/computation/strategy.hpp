#ifndef COMPUTATION_STRATEGY_HPP
#define COMPUTATION_STRATEGY_HPP

#include "grid/partition.hpp"
#include "util/cell.hpp"
#include "parameters.hpp"

namespace computation {

class strategy {
    public:
        strategy(index_grid_type const& index_grid, parameters const& params) : index(index_grid), p(params) {}

        virtual void set_velocity_on_boundary(vector_grid_type& uv_grid) {}
        virtual void set_pressure_on_boundary(scalar_grid_type& p_grid) {}
        virtual void compute_fg(vector_grid_type& fg_grid, vector_grid_type const& uv_grid, RealType dt) {}
        virtual void compute_rhs(scalar_grid_type& rhs_grid, vector_grid_type const& fg_grid, RealType dt) {}

        virtual void sor_cycle(scalar_grid_type& p_grid, scalar_grid_type const& rhs_grid) {}
        virtual hpx::future<RealType> compute_residual(scalar_grid_type const& p_grid, scalar_grid_type const& rhs_grid) {}

        virtual hpx::future<std::pair<RealType, RealType> > update_velocities(vector_grid_type& uv_grid,
                                                        vector_grid_type const& fg_grid, scalar_grid_type const& p_grid, RealType dt) {}

    protected:
        index_grid_type const& index;
        parameters const& p;

        inline uint get_index(uint k, uint l) { return l * p.num_partitions_x + k;}

};

}//namespace

#endif
