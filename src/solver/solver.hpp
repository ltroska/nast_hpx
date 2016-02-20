#ifndef SOLVER_SOLVER_HPP
#define SOLVER_SOLVER_HPP

#include "grid/partition.hpp"
#include "parameters.hpp"

namespace solver {

class solver {
    protected:
        typedef std::vector<grid::partition> grid_type;
        typedef std::vector<std::pair<uint, uint> > index_grid_type;

    public:
        solver(index_grid_type const& index_grid, parameters const& params) : index(index_grid), p(params) {}
        virtual void set_velocity_on_boundary(grid_type& u_grid, grid_type& v_grid) {}

    protected:
        index_grid_type const& index;
        parameters const& p;

};

}//namespace

#endif
