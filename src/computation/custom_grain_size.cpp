#include "custom_grain_size.hpp"

#include "util/helpers.hpp"

namespace computation {

vector_partition set_velocity_on_boundary_action(hpx::shared_future<vector_data> center_fut, hpx::naming::id_type const where,
                                                    uint global_i, uint global_j, uint i_max, uint j_max)
{
    vector_data center = center_fut.get();

    uint size_x = center.size_x();
    uint size_y = center.size_y();

    bool is_left = (global_i == 0);
    bool is_right = (global_i + size_x > i_max);

    bool is_bottom = (global_j == 0);
    bool is_top = (global_j + size_y > j_max);

    uint start_i = (is_left ? 1 : 0);
    uint end_i = (is_right ? size_x - 1 : size_x);
    uint start_j = (is_bottom ? 1 : 0);
    uint end_j = (is_top ? size_y - 1 : size_y);

    //eq 17+18
    if (is_left)
    {
        for (uint j = start_j; j < end_j; j++)
        {
            vector_cell& cell = center.get_cell_ref(0, j);
            vector_cell cell2 = center.get_cell(1, j);

            cell.first = 0;
            cell.second = -cell2.second;
        }
    }

    if (is_right)
    {
        for (uint j = start_j; j < end_j; j++)
        {
            vector_cell& cell = center.get_cell_ref(size_x - 2, j);
            vector_cell& cell2 = center.get_cell_ref(size_x - 1, j);

            cell.first = 0;
            cell2.second = -cell.second;
        }
    }

    if (is_bottom)
    {
        for (uint i = start_i; i < end_i; i++)
        {
            vector_cell& cell = center.get_cell_ref(i, 0);
            vector_cell cell2 = center.get_cell(i, 1);

            cell.second = 0;
            cell.first = -cell2.first;

        }
    }

    if (is_top)
    {
        for (uint i = start_i; i < end_i; i++)
        {
            vector_cell& cell = center.get_cell_ref(i, size_y - 2);
            vector_cell& cell2 = center.get_cell_ref(i, size_y - 1);

            cell2.second = 0;
            cell2.first = 2. - cell.first;
        }
    }

    return vector_partition(where, center);
}

vector_partition dispatch_set_velocity_on_boundary(vector_partition const& center, uint global_i, uint global_j, uint i_max, uint j_max)
{
    hpx::shared_future<vector_data> center_data = center.get_data(CENTER);

    hpx::naming::id_type const where = center.get_id();

    return hpx::dataflow(
            hpx::launch::async,
            &set_velocity_on_boundary_action,
            center_data,
            where,
            global_i,
            global_j,
            i_max,
            j_max
    );
}

void custom_grain_size::set_velocity_on_boundary(vector_grid_type& uv_grid)
{
    for (uint l = 1; l < p.num_partitions_y - 1; l++)
    {
        for (uint k = 1; k < p.num_partitions_x - 1; k++)
        {
            //skip interior
            if (!(k == 1 || k == p.num_partitions_x - 2 || l == 1 || l == p.num_partitions_y - 2))
                continue;

            vector_partition& next = uv_grid[get_index(k, l)];

            next =
                hpx::dataflow(
                    hpx::launch::async,
                    &dispatch_set_velocity_on_boundary,
                    next,
                    index[get_index(k, l)].first,
                    index[get_index(k, l)].second,
                    p.i_max,
                    p.j_max
            );
        }
    }

}

}//namespace

