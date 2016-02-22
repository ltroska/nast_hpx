#include "custom_grain_size.hpp"

#include "util/helpers.hpp"

namespace computation {

vector_partition set_velocity_on_boundary_action(hpx::naming::id_type const where, hpx::shared_future<vector_data> center_fut,
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
            where,
            center_data,
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

vector_partition compute_fg_action(hpx::naming::id_type const where, hpx::shared_future<vector_data> center,
                                        hpx::shared_future<vector_data> left, hpx::shared_future<vector_data> right,
                                        hpx::shared_future<vector_data> bottom, hpx::shared_future<vector_data> top,
                                        hpx::shared_future<vector_data> bottomright, hpx::shared_future<vector_data> topleft,
                                        uint global_i, uint global_j, uint i_max, uint j_max, RealType dx, RealType dy,
                                        RealType re, RealType alpha, RealType dt)
{

    vector_data uv_center = center.get();
    vector_data uv_left = left.get();
    vector_data uv_right = right.get();
    vector_data uv_bottom = bottom.get();
    vector_data uv_top = top.get();
    vector_data uv_bottomright = bottomright.get();
    vector_data uv_topleft = topleft.get();

    return vector_partition(where, uv_center);
}

vector_partition dispatch_compute_fg(vector_partition const& uv_center, vector_partition const& uv_left, vector_partition const& uv_right,
                                        vector_partition const& uv_bottom, vector_partition const& uv_top, vector_partition const& uv_bottomright,
                                        vector_partition const& uv_topleft, uint global_i, uint global_j, uint i_max, uint j_max, RealType dx,
                                        RealType dy, RealType re, RealType alpha, RealType dt)
{
    hpx::shared_future<vector_data> uv_center_data = uv_center.get_data(CENTER);
    hpx::shared_future<vector_data> uv_left_data = uv_left.get_data(LEFT);
    hpx::shared_future<vector_data> uv_right_data = uv_right.get_data(RIGHT);
    hpx::shared_future<vector_data> uv_bottom_data = uv_bottom.get_data(BOTTOM);
    hpx::shared_future<vector_data> uv_top_data = uv_top.get_data(TOP);
    hpx::shared_future<vector_data> uv_bottomright_data = uv_bottomright.get_data(BOTTOM_RIGHT);
    hpx::shared_future<vector_data> uv_topleft_data = uv_topleft.get_data(TOP_LEFT);

    hpx::naming::id_type const where = uv_center.get_id();

    return hpx::dataflow(
            hpx::launch::async,
            &compute_fg_action,
            where,
            uv_center_data, uv_left_data, uv_right_data,
            uv_bottom_data, uv_top_data, uv_bottomright_data, uv_topleft_data,
            global_i,
            global_j,
            i_max, j_max, dx, dy, re, alpha, dt
    );

}

void custom_grain_size::compute_fg(vector_grid_type& fg_grid, vector_grid_type const& uv_grid, RealType dt)
{
    for (uint l = 1; l < p.num_partitions_y - 1; l++)
    {
        for (uint k = 1; k < p.num_partitions_x - 1; k++)
        {

            vector_partition& next = fg_grid[get_index(k, l)];

            next =
                hpx::dataflow(
                    hpx::launch::async,
                    &dispatch_compute_fg,
                    uv_grid[get_index(k, l)], //center
                    uv_grid[get_index(k-1, l)], //left
                    uv_grid[get_index(k+1, l)], //right
                    uv_grid[get_index(k, l-1)], //bottom
                    uv_grid[get_index(k, l+1)], //top
                    uv_grid[get_index(k+1, l-1)], //bottomright
                    uv_grid[get_index(k-1, l+1)], //topleft
                    index[get_index(k, l)].first,
                    index[get_index(k, l)].second,
                    p.i_max, p.j_max, p.dx, p.dy, p.re, p.alpha, dt
            );
        }
    }
}

}//namespace

