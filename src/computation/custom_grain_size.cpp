#include "custom_grain_size.hpp"

#include "util/helpers.hpp"
#include "stencils.hpp"

namespace computation {

/*
// ------------------------------------------------------------ SET VELOCITY ON BOUNDARY ------------------------------------------------------------ //
*/

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

/*
// ------------------------------------------------------------ COMPUTE FG ------------------------------------------------------------ //
*/

vector_partition compute_fg_action(hpx::naming::id_type const where, hpx::shared_future<vector_data> center_fut,
                                        hpx::shared_future<vector_data> left_fut, hpx::shared_future<vector_data> right_fut,
                                        hpx::shared_future<vector_data> bottom_fut, hpx::shared_future<vector_data> top_fut,
                                        hpx::shared_future<vector_data> bottomright_fut, hpx::shared_future<vector_data> topleft_fut,
                                        uint global_i, uint global_j, uint i_max, uint j_max, RealType dx, RealType dy,
                                        RealType re, RealType alpha, RealType dt)
{

    vector_data uv_center = center_fut.get();
    vector_data uv_left = left_fut.get();
    vector_data uv_right = right_fut.get();
    vector_data uv_bottom = bottom_fut.get();
    vector_data uv_top = top_fut.get();
    vector_data uv_bottomright = bottomright_fut.get();
    vector_data uv_topleft = topleft_fut.get();

    uint size_x = uv_center.size_x();
    uint size_y = uv_center.size_y();

    vector_data fg_data(size_x, size_y);

    bool is_left = (global_i == 0);
    bool is_right = (global_i + size_x > i_max);

    bool is_bottom = (global_j == 0);
    bool is_top = (global_j + size_y > j_max);

    //do computation for i = 1, ..., i_max-1, j = 1, ..., j_max-1 first
    uint start_i = (is_left ? 1 : 0);
    uint end_i = (is_right ? size_x - 2 : size_x);
    uint start_j = (is_bottom ? 1 : 0);
    uint end_j = (is_top ? size_y - 2 : size_y);

    for (uint j = start_j; j < end_j; j++)
        for (uint i = start_i; i < end_i; i++)
        {
            vector_cell& fg_cell = fg_data.get_cell_ref(i, j);

            vector_cell const center = uv_center.get_cell(i, j);
            vector_cell const left = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, LEFT);
            vector_cell const right = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, RIGHT);
            vector_cell const bottom = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, BOTTOM);
            vector_cell const top = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, TOP);
            vector_cell const bottomright = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, BOTTOM_RIGHT);
            vector_cell const topleft = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, TOP_LEFT);

            fg_cell.first = center.first + dt * (
                                    1./re * (second_derivative_fwd_bkwd_x(right.first, center.first, left.first, dx)
                                                + second_derivative_fwd_bkwd_y(top.first, center.first, bottom.first, dy))

                                    - first_derivative_of_square_x(right.first, center.first, left.first, dx, alpha)
                                    - first_derivative_of_product_y(right.second, center.second, bottom.second, bottomright.second,
                                                                    bottom.first, center.first, top.first, dy, alpha)
                                    );

            fg_cell.second = center.second + dt * (
                                    1./re * (second_derivative_fwd_bkwd_x(right.second, center.second, left.second, dx)
                                                + second_derivative_fwd_bkwd_y(top.second, center.second, bottom.second, dy))

                                    - first_derivative_of_product_x(left.first, center.first, top.first, topleft.first,
                                                                    left.second, center.second, right.second, dx, alpha)
                                    - first_derivative_of_square_y(top.second, center.second, bottom.second, dy, alpha)
                                    );
            }

    // compute top strip i = 1, ..., i_max-1, j = j_max
    if (is_top)
    {
        uint j = size_y - 2;

        for (uint i = start_i; i < end_i; i++)
        {
            vector_cell& fg_cell = fg_data.get_cell_ref(i, j);

            vector_cell const center = uv_center.get_cell(i, j);
            vector_cell const left = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, LEFT);
            vector_cell const right = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, RIGHT);
            vector_cell const bottom = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, BOTTOM);
            vector_cell const top = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, TOP);
            vector_cell const bottomright = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, BOTTOM_RIGHT);
            vector_cell const topleft = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, TOP_LEFT);

            fg_cell.first = center.first + dt * (
                                    1./re * (second_derivative_fwd_bkwd_x(right.first, center.first, left.first, dx)
                                                + second_derivative_fwd_bkwd_y(top.first, center.first, bottom.first, dy))

                                    - first_derivative_of_square_x(right.first, center.first, left.first, dx, alpha)
                                    - first_derivative_of_product_y(right.second, center.second, bottom.second, bottomright.second,
                                                                    bottom.first, center.first, top.first, dy, alpha)
                                    );
        }
    }

    // compute right strip i = i_max, j = 1, ..., j_max-1
    if (is_right)
    {
        uint i = size_x - 2;

        for (uint j = start_j; j < end_j; j++)
        {
            vector_cell& fg_cell = fg_data.get_cell_ref(i, j);

            vector_cell const center = uv_center.get_cell(i, j);
            vector_cell const left = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, LEFT);
            vector_cell const right = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, RIGHT);
            vector_cell const bottom = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, BOTTOM);
            vector_cell const top = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, TOP);
            vector_cell const bottomright = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, BOTTOM_RIGHT);
            vector_cell const topleft = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, TOP_LEFT);

            fg_cell.second = center.second + dt * (
                                    1./re * (second_derivative_fwd_bkwd_x(right.second, center.second, left.second, dx)
                                                + second_derivative_fwd_bkwd_y(top.second, center.second, bottom.second, dy))

                                    - first_derivative_of_product_x(left.first, center.first, top.first, topleft.first,
                                                                    left.second, center.second, right.second, dx, alpha)
                                    - first_derivative_of_square_y(top.second, center.second, bottom.second, dy, alpha)
                                    );
        }
    }

    return vector_partition(where, fg_data);
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

/*
// ------------------------------------------------------------ COMPUTE RHS ------------------------------------------------------------ //
*/
scalar_partition compute_rhs_action(hpx::naming::id_type const where, hpx::shared_future<vector_data> center_fut,
                                        hpx::shared_future<vector_data> left_fut,  hpx::shared_future<vector_data> bottom_fut,
                                        uint global_i, uint global_j, uint i_max, uint j_max, RealType dx, RealType dy, RealType dt)
{

    vector_data fg_center = center_fut.get();
    vector_data fg_left = left_fut.get();
    vector_data fg_bottom = bottom_fut.get();

    uint size_x = fg_center.size_x();
    uint size_y = fg_center.size_y();

    scalar_data rhs_data(size_x, size_y);

    bool is_left = (global_i == 0);
    bool is_right = (global_i + size_x > i_max);

    bool is_bottom = (global_j == 0);
    bool is_top = (global_j + size_y > j_max);

    //do computation for i = 1, ..., i_max, j = 1, ..., j_max
    uint start_i = (is_left ? 1 : 0);
    uint end_i = (is_right ? size_x - 1 : size_x);
    uint start_j = (is_bottom ? 1 : 0);
    uint end_j = (is_top ? size_y - 1 : size_y);


    for (uint j = start_j; j < end_j; j++)
        for (uint i = start_i; i < end_i; i++)
        {
            vector_cell const center = fg_center.get_cell(i, j);
            vector_cell const left = get_neighbor_cell(fg_center, fg_left, fg_left, fg_bottom, fg_bottom, fg_bottom, fg_bottom, fg_bottom, fg_bottom, i, j, LEFT);
            vector_cell const bottom = get_neighbor_cell(fg_center, fg_left, fg_left, fg_bottom, fg_bottom, fg_bottom, fg_bottom, fg_bottom, fg_bottom, i, j, BOTTOM);

            rhs_data.get_cell_ref(i, j).value = 1./dt * ( (center.first - left.first)/dx + (center.second - bottom.second)/dy);
        }

    return scalar_partition(where, rhs_data);
}

scalar_partition dispatch_compute_rhs(vector_partition const& fg_center, vector_partition const& fg_left, vector_partition const& fg_bottom,
                                        uint global_i, uint global_j, uint i_max, uint j_max, RealType dx, RealType dy, RealType dt)
{
    hpx::shared_future<vector_data> fg_center_data = fg_center.get_data(CENTER);
    hpx::shared_future<vector_data> fg_left_data = fg_left.get_data(LEFT);
    hpx::shared_future<vector_data> fg_bottom_data = fg_bottom.get_data(BOTTOM);

    hpx::naming::id_type const where = fg_center.get_id();

    return hpx::dataflow(
            hpx::launch::async,
            &compute_rhs_action,
            where,
            fg_center_data, fg_left_data, fg_bottom_data,
            global_i,
            global_j,
            i_max, j_max, dx, dy, dt
    );

}

void custom_grain_size::compute_rhs(scalar_grid_type& rhs_grid, vector_grid_type const& fg_grid, RealType dt)
{
    for (uint l = 1; l < p.num_partitions_y - 1; l++)
    {
        for (uint k = 1; k < p.num_partitions_x - 1; k++)
        {

            scalar_partition& next = rhs_grid[get_index(k, l)];

            next =
                hpx::dataflow(
                    hpx::launch::async,
                    &dispatch_compute_rhs,
                    fg_grid[get_index(k, l)], //center
                    fg_grid[get_index(k-1, l)], //left
                    fg_grid[get_index(k, l-1)], //bottom
                    index[get_index(k, l)].first,
                    index[get_index(k, l)].second,
                    p.i_max, p.j_max, p.dx, p.dy, dt
            );
        }
    }
}


}//namespace

