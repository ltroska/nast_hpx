#include "custom_chunk_solver.hpp"

#include "util/helpers.hpp"

namespace solver {

template<typename Func, typename... Args>
auto dispatch(Func&& f, grid::partition const& p, Args&&... args)
-> decltype(f(std::declval<hpx::future<grid::partition_data<> > >(), std::declval<hpx::naming::id_type>(), args...))
{
   return hpx::dataflow(
            hpx::launch::async,
            std::forward<Func>(f),
            p.get_data(CENTER),
            p.get_id(),
            std::forward<Args>(args)...
    );
}


grid::partition set_u_on_boundary(hpx::future<grid::partition_data<> > p, const hpx::naming::id_type where, uint global_i, uint global_j, uint i_max, uint j_max)
{
    grid::partition_data<> current = p.get();

    uint size_x = current.size_x();
    uint size_y = current.size_y();

    for (uint j = 0; j < size_y; j++)
    {
        if (in_range(0, 0, 1, j_max, global_i, global_j + j))
            current.get_cell_ref(0, j) = 1;
    }

    return grid::partition(where, current);
}

grid::partition set_v_on_boundary(hpx::future<grid::partition_data<> > p, const hpx::naming::id_type where, uint global_i, uint global_j, uint i_max, uint j_max)
{
    grid::partition_data<> current = p.get();

    uint size_x = current.size_x();
    uint size_y = current.size_y();

    for (uint j = 0; j < size_y; j++)
    {
        if (in_range(0, 0, 1, j_max, global_i, global_j + j))
            current.get_cell_ref(0, j) = 6;
    }

    return grid::partition(where, current);
}

void custom_chunk_solver::set_velocity_on_boundary(grid_type& u_grid, grid_type& v_grid)
{
/*
*@TODO: maybe add bool is_*_bnd !
*/
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
    END_FOR
}

}//namespace
