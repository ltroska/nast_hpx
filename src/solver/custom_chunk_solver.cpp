#include "custom_chunk_solver.hpp"

#include "util/helpers.hpp"

namespace solver {

void custom_chunk_solver::set_velocity_on_boundary(vector_grid_type& uv_grid)
{
/*
*@TODO: maybe add bool is_*_bnd !
*/
/*
    FOR_EVERY_BOUNDARY_PARTITION
    {
            u_grid[p.num_partitions_x * l + k] = dispatch(
                                                    set_u_on_boundary,
                                                    u_grid[p.num_partitions_x * l + k],
                                                    index[l * p.num_partitions_x + k].first,
                                                    index[l * p.num_partitions_x + k].second,
                                                    p.i_max,
                                                    p.j_max
                                                    );

             v_grid[p.num_partitions_x * l + k] = dispatch(
                                                    set_v_on_boundary,
                                                    v_grid[p.num_partitions_x * l + k],
                                                    index[l * p.num_partitions_x + k].first,
                                                    index[l * p.num_partitions_x + k].second,
                                                    p.i_max,
                                                    p.j_max
                                                    );
    }
    END_FOR*/
}

}//namespace
