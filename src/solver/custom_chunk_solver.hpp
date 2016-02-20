#ifndef SOLVER_CUSTOM_CHUNK_SOLVER_HPP
#define SOLVER_CUSTOM_CHUNK_SOLVER_HPP

#include "solver.hpp"

namespace solver {

class custom_chunk_solver : public solver {

    using solver::solver;

    virtual void set_velocity_on_boundary(vector_grid_type& uv_grid);
};

}//namespace

#endif
