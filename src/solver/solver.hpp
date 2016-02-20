#ifndef SOLVER_SOLVER_HPP
#define SOLVER_SOLVER_HPP

#include "grid/partition.hpp"
#include "util/cell.hpp"
#include "parameters.hpp"

namespace solver {

class solver {
    protected:
        typedef grid::partition<scalar_cell> scalar_partition;
        typedef grid::partition<vector_cell> vector_partition;
        typedef std::vector<scalar_partition> scalar_grid_type;
        typedef std::vector<vector_partition> vector_grid_type;
        typedef std::vector<std::pair<uint, uint> > index_grid_type;

    public:
        solver(index_grid_type const& index_grid, parameters const& params) : index(index_grid), p(params) {}
        virtual void set_velocity_on_boundary(vector_grid_type& uv_grid) {}

    protected:
        index_grid_type const& index;
        parameters const& p;

};

}//namespace

#endif
