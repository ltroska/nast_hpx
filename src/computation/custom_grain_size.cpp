#include "custom_grain_size.hpp"

#include "util/helpers.hpp"
#include "stencils.hpp"
#include <hpx/parallel/algorithm.hpp>

namespace computation {

/*
// ------------------------------------------------------------ SET VELOCITY ON BOUNDARY ------------------------------------------------------------ //
*/

vector_partition set_velocity_on_boundary_action(hpx::naming::id_type const where, hpx::shared_future<vector_data> center_fut,
                                                    uint global_i, uint global_j, uint i_max, uint j_max)
{
    /*
    *@TODO: maybe create new vector_data here
    */
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
            vector_cell const cell2 = center.get_cell(1, j);

            cell.first = 0;
            cell.second = -cell2.second;
        }
    }

    if (is_bottom)
    {
        for (uint i = start_i; i < end_i; i++)
        {
            vector_cell& cell = center.get_cell_ref(i, 0);
            vector_cell const cell2 = center.get_cell(i, 1);

            cell.second = 0;
            cell.first = -cell2.first;

        }
    }

    // need to buffer this value, in case this is top and(!) right boundary, because then
    // the right and the top equations both try set this value, resulting in an inconsistent state
    RealType conflicting_u = center.get_cell(size_x - 2, size_y - 2).first;

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

    if (is_top)
    {
        for (uint i = start_i; i < end_i; i++)
        {
            vector_cell& cell = center.get_cell_ref(i, size_y - 2);
            vector_cell& cell2 = center.get_cell_ref(i, size_y - 1);

            cell.second = 0;
            cell2.first = 2. - cell.first;
        }
    }

    if (is_right && is_top)
        center.get_cell_ref(size_x - 2, size_y - 1).first = 2 - conflicting_u;

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

    // compute top strip i = 1, ..., i_max-1, j = j_max for F and set top boundary G = v
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

            fg_cell.second = center.second;
        }
    }

    // compute right strip i = i_max, j = 1, ..., j_max-1 and set right boundary F = u
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

            fg_cell.first = center.first;
        }
    }

    end_i = (is_right ? size_x - 1 : size_x);
    end_j = (is_top ? size_y - 1 : size_y);

    // setting left boundary F = u
    if (is_left)
    {
        uint i = 0;
        for (uint j = start_j; j < end_j; j++)
        {
            vector_cell& fg_cell = fg_data.get_cell_ref(i, j);
            vector_cell const center = uv_center.get_cell(i, j);

            fg_cell.first = center.first;
        }
    }

    // setting bottom boundary G = v
    if (is_bottom)
    {
        uint j = 0;
        for (uint i = start_i; i < end_i; i++)
        {
            vector_cell& fg_cell = fg_data.get_cell_ref(i, j);
            vector_cell const center = uv_center.get_cell(i, j);

            fg_cell.second = center.second;
        }
    }

    // set cells we missed in above computation
    if (is_top)
        fg_data.get_cell_ref(size_x - 2, size_y - 2).second = uv_center.get_cell(size_x - 2, size_y - 2).second;

    if (is_right)
        fg_data.get_cell_ref(size_x - 2, size_y - 2).first = uv_center.get_cell(size_x - 2, size_y - 2).first;

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

/*
// ------------------------------------------------------------ SET VELOCITY ON BOUNDARY ------------------------------------------------------------ //
*/

scalar_partition set_pressure_on_boundary_action(hpx::naming::id_type const where, hpx::shared_future<scalar_data> center_fut,
                                                    uint global_i, uint global_j, uint i_max, uint j_max)
{
    scalar_data center = center_fut.get();

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

    if (is_left)
    {
        for (uint j = start_j; j < end_j; j++)
        {
            scalar_cell& cell = center.get_cell_ref(0, j);
            scalar_cell const cell2 = center.get_cell(1, j);

            cell.value = cell2.value;
        }
    }

    if (is_right)
    {
        for (uint j = start_j; j < end_j; j++)
        {
            scalar_cell& cell = center.get_cell_ref(size_x - 1, j);
            scalar_cell const cell2 = center.get_cell_ref(size_x - 2, j);

            cell.value = cell2.value;
        }
    }

    if (is_bottom)
    {
        for (uint i = start_i; i < end_i; i++)
        {
            scalar_cell& cell = center.get_cell_ref(i, 0);
            scalar_cell const cell2 = center.get_cell(i, 1);

            cell.value = cell2.value;
        }
    }

    if (is_top)
    {
        for (uint i = start_i; i < end_i; i++)
        {
            scalar_cell& cell = center.get_cell_ref(i, size_y - 1);
            scalar_cell const cell2 = center.get_cell_ref(i, size_y - 2);

            cell.value = cell2.value;
        }
    }

    return scalar_partition(where, center);
}

scalar_partition dispatch_set_pressure_on_boundary(scalar_partition const& center, uint global_i, uint global_j, uint i_max, uint j_max)
{
    hpx::shared_future<scalar_data> center_data = center.get_data(CENTER);

    hpx::naming::id_type const where = center.get_id();

    return hpx::dataflow(
            hpx::launch::async,
            &set_pressure_on_boundary_action,
            where,
            center_data,
            global_i,
            global_j,
            i_max,
            j_max
    );
}

void custom_grain_size::set_pressure_on_boundary(scalar_grid_type& p_grid)
{
    for (uint l = 1; l < p.num_partitions_y - 1; l++)
    {
        for (uint k = 1; k < p.num_partitions_x - 1; k++)
        {
            //skip interior
            if (!(k == 1 || k == p.num_partitions_x - 2 || l == 1 || l == p.num_partitions_y - 2))
                continue;

            scalar_partition& next = p_grid[get_index(k, l)];

            next =
                hpx::dataflow(
                    hpx::launch::async,
                    &dispatch_set_pressure_on_boundary,
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
// ------------------------------------------------------------ SOR CYCLE ------------------------------------------------------------ //
*/

scalar_partition sor_cycle_action(hpx::naming::id_type const where, hpx::shared_future<scalar_data> center_fut, hpx::shared_future<scalar_data> left_fut,
                                        hpx::shared_future<scalar_data> right_fut, hpx::shared_future<scalar_data> bottom_fut, hpx::shared_future<scalar_data> top_fut,
                                        hpx::shared_future<scalar_data> rhs_fut, uint global_i, uint global_j, uint i_max, uint j_max,
                                        RealType omega, RealType dx, RealType dy)
{
    /*
    *@TODO: maybe create new scalar_data
    */

    scalar_data p_center = center_fut.get();
    scalar_data p_left = left_fut.get();
    scalar_data p_right = right_fut.get();
    scalar_data p_bottom = bottom_fut.get();
    scalar_data p_top = top_fut.get();
    scalar_data rhs_center = rhs_fut.get();

    uint size_x = p_center.size_x();
    uint size_y = p_center.size_y();

    scalar_data center(p_center);

    bool is_left = (global_i == 0);
    bool is_right = (global_i + size_x > i_max);
    bool is_bottom = (global_j == 0);
    bool is_top = (global_j + size_y > j_max);

    uint start_i = (is_left ? 1 : 0);
    uint end_i = (is_right ? size_x - 1 : size_x);
    uint start_j = (is_bottom ? 1 : 0);
    uint end_j = (is_top ? size_y - 1 : size_y);


    RealType over_dx_sq = 1./std::pow(dx, 2);
    RealType over_dy_sq = 1./std::pow(dy, 2);
    RealType part1 = 1. - omega;
    RealType part2 = omega / (2. * (over_dx_sq + over_dy_sq));

  /*  auto range = boost::irange(0, (int)(size_x*size_y));

    hpx::parallel::for_each(hpx::parallel::par, boost::begin(range),
        boost::end(range),
        [&](uint cnt)
        {
            uint i = cnt%size_x;
            uint j = cnt/size_x;
            scalar_cell& next_p = p_center.get_cell_ref(i, j);
            scalar_cell const current_rhs = rhs_center.get_cell(i, j);
            scalar_cell const left = get_left_neighbor(p_center, p_left, i, j);
            scalar_cell const right = get_right_neighbor(p_center, p_right, i, j);
            scalar_cell const bottom = get_bottom_neighbor(p_center, p_bottom, i, j);
            scalar_cell const top = get_top_neighbor(p_center, p_top, i, j);

            next_p.value = part1 * next_p.value
                            + part2 * ( (right.value + left.value)*over_dx_sq + (top.value + bottom.value)*over_dy_sq - current_rhs.value);
        }
    );*/

    for (uint j = start_j; j < end_j; j ++)
    {
        for (uint i = start_i; i < end_i; i++)
        {
            scalar_cell& next_p = center.get_cell_ref(i, j);
            scalar_cell const current_rhs = rhs_center.get_cell(i, j);
            scalar_cell const left = get_left_neighbor(center, p_left, i, j);
            scalar_cell const right = get_right_neighbor(center, p_right, i, j);
            scalar_cell const bottom = get_bottom_neighbor(center, p_bottom, i, j);
            scalar_cell const top = get_top_neighbor(center, p_top, i, j);

            next_p.value = part1 * next_p.value
                            + part2 * ( (right.value + left.value)*over_dx_sq + (top.value + bottom.value)*over_dy_sq - current_rhs.value);
        }
    }

    return scalar_partition(where, center);
}

scalar_partition dispatch_sor_cycle(scalar_partition const& center, scalar_partition const& left, scalar_partition const& right,
                                        scalar_partition const& bottom, scalar_partition const& top, scalar_partition const& rhs,
                                        uint global_i, uint global_j, uint i_max, uint j_max, RealType omega, RealType dx, RealType dy)
{
    hpx::shared_future<scalar_data> center_data = center.get_data(CENTER);
    hpx::shared_future<scalar_data> left_data = left.get_data(LEFT);
    hpx::shared_future<scalar_data> right_data = right.get_data(RIGHT);
    hpx::shared_future<scalar_data> bottom_data = bottom.get_data(BOTTOM);
    hpx::shared_future<scalar_data> top_data = top.get_data(TOP);
    hpx::shared_future<scalar_data> rhs_data = rhs.get_data(CENTER);

    hpx::naming::id_type const where = center.get_id();

    return hpx::dataflow(
            hpx::launch::async,
            &sor_cycle_action,
            where,
            center_data, left_data, right_data,
            bottom_data, top_data, rhs_data,
            global_i,
            global_j,
            i_max,
            j_max,
            omega,
            dx,
            dy
    );
}


RealType compute_residual_action(hpx::shared_future<scalar_data> center_fut, hpx::shared_future<scalar_data> left_fut,
                            hpx::shared_future<scalar_data> right_fut, hpx::shared_future<scalar_data> bottom_fut, hpx::shared_future<scalar_data> top_fut,
                            hpx::shared_future<scalar_data> rhs_fut, uint global_i, uint global_j, uint i_max, uint j_max, RealType dx, RealType dy)
{
    scalar_data p_center = center_fut.get();
    scalar_data p_left = left_fut.get();
    scalar_data p_right = right_fut.get();
    scalar_data p_bottom = bottom_fut.get();
    scalar_data p_top = top_fut.get();
    scalar_data rhs_center = rhs_fut.get();

    uint size_x = p_center.size_x();
    uint size_y = p_center.size_y();

    bool is_left = (global_i == 0);
    bool is_right = (global_i + size_x > i_max);
    bool is_bottom = (global_j == 0);
    bool is_top = (global_j + size_y > j_max);

    uint start_i = (is_left ? 1 : 0);
    uint end_i = (is_right ? size_x - 1 : size_x);
    uint start_j = (is_bottom ? 1 : 0);
    uint end_j = (is_top ? size_y - 1 : size_y);

    RealType local_residual = 0;

    RealType over_dx_sq = 1./(dx*dx);
    RealType over_dy_sq = 1./(dy*dy);

    for (uint j = start_j; j < end_j; j++)
        for (uint i = start_i; i < end_i; i++)
        {
            scalar_cell const center = p_center.get_cell(i, j);
            scalar_cell const rhs = rhs_center.get_cell(i, j);
            scalar_cell const left = get_neighbor_cell(p_center, p_left, p_right, p_bottom, p_top, p_top, p_top, p_top, p_top, i, j, LEFT);
            scalar_cell const right = get_neighbor_cell(p_center, p_left, p_right, p_bottom, p_top, p_top, p_top, p_top, p_top, i, j, RIGHT);
            scalar_cell const bottom = get_neighbor_cell(p_center, p_left, p_right, p_bottom, p_top, p_top, p_top, p_top, p_top, i, j, BOTTOM);
            scalar_cell const top = get_neighbor_cell(p_center, p_left, p_right, p_bottom, p_top, p_top, p_top, p_top, p_top, i, j, TOP);

            RealType tmp = (right.value - 2*center.value + left.value)*over_dx_sq + (top.value - 2*center.value + bottom.value)*over_dy_sq - rhs.value;
            local_residual += std::pow(tmp, 2);
        }

    return local_residual/(i_max*j_max);
}

hpx::future<RealType> dispatch_compute_residual(scalar_partition const& center, scalar_partition const& left, scalar_partition const& right,
                                        scalar_partition const& bottom, scalar_partition const& top, scalar_partition const& rhs,
                                        uint global_i, uint global_j, uint i_max, uint j_max, RealType dx, RealType dy)
{
    hpx::shared_future<scalar_data> center_data = center.get_data(CENTER);
    hpx::shared_future<scalar_data> left_data = left.get_data(LEFT);
    hpx::shared_future<scalar_data> right_data = right.get_data(RIGHT);
    hpx::shared_future<scalar_data> bottom_data = bottom.get_data(BOTTOM);
    hpx::shared_future<scalar_data> top_data = top.get_data(TOP);
    hpx::shared_future<scalar_data> rhs_data = rhs.get_data(CENTER);


    return hpx::dataflow(
            hpx::launch::async,
            &compute_residual_action,
            center_data, left_data, right_data,
            bottom_data, top_data, rhs_data,
            global_i,
            global_j,
            i_max,
            j_max,
            dx,
            dy
    );
}

hpx::future<RealType> custom_grain_size::sor_cycle(scalar_grid_type& p_grid, scalar_grid_type const& rhs_grid)
{
    // odd cells first
    for (uint l = 1; l < p.num_partitions_y - 1; l++)
    {
        for (uint k = 1; k < p.num_partitions_x - 1; k++)
        {

            p_grid[get_index(k, l)] =
                hpx::dataflow(
                    hpx::launch::async,
                    &dispatch_sor_cycle,
                    p_grid[get_index(k, l)],
                    p_grid[get_index(k-1, l)], //left
                    p_grid[get_index(k+1, l)], //right
                    p_grid[get_index(k, l-1)], //bottom
                    p_grid[get_index(k, l+1)], //top
                    rhs_grid[get_index(k, l)],
                    index[get_index(k, l)].first,
                    index[get_index(k, l)].second,
                    p.i_max,
                    p.j_max,
                    p.omega,
                    p.dx,
                    p.dy
            );
        }
    }

    // residuals
    hpx::future<RealType> residual = hpx::make_ready_future(0.0);
    for (uint l = 1; l < p.num_partitions_y - 1; l++)
    {
        for (uint k = 1; k < p.num_partitions_x - 1; k++)
        {
                residual = hpx::dataflow(
                    hpx::launch::async,
                    [](hpx::future<RealType> prev_sum, hpx::future<RealType> next_summand)
                        -> RealType
                    {
                        return prev_sum.get() + next_summand.get();
                    }
                    , residual
                    , hpx::dataflow(
                            hpx::launch::async,
                            &dispatch_compute_residual,
                            p_grid[get_index(k, l)], //center
                            p_grid[get_index(k-1, l)], //left
                            p_grid[get_index(k+1, l)], //right
                            p_grid[get_index(k, l-1)], //bottom
                            p_grid[get_index(k, l+1)], //top
                            rhs_grid[get_index(k, l)],
                            index[get_index(k, l)].first,
                            index[get_index(k, l)].second,
                            p.i_max,
                            p.j_max,
                            p.dx,
                            p.dy
                        )
                    );
        }
    }

    return residual;
}

/*
// ------------------------------------------------------------ UPDATE VELOCITIES ------------------------------------------------------------ //
*/

vector_partition update_velocities_action(hpx::naming::id_type const where, hpx::shared_future<vector_data> uv_center_fut, hpx::shared_future<scalar_data> center_fut,
                            hpx::shared_future<scalar_data> right_fut, hpx::shared_future<scalar_data> top_fut, hpx::shared_future<vector_data> fg_fut,
                            uint global_i, uint global_j, uint i_max, uint j_max, RealType dx, RealType dy, RealType dt)
{
    /*
    *@TODO: Maybe create new vector_data here
    */

    vector_data uv_center = uv_center_fut.get();
    scalar_data p_center = center_fut.get();
    scalar_data p_right = right_fut.get();
    scalar_data p_top = top_fut.get();
    vector_data fg_center = fg_fut.get();

    uint size_x = uv_center.size_x();
    uint size_y = uv_center.size_y();

    bool is_left = (global_i == 0);
    bool is_right = (global_i + size_x > i_max);
    bool is_bottom = (global_j == 0);
    bool is_top = (global_j + size_y > j_max);

    RealType over_dx = 1./dx;
    RealType over_dy = 1./dy;

    //do computation for i = 1, ..., i_max-1, j = 1, ..., j_max-1 first
    uint start_i = (is_left ? 1 : 0);
    uint end_i = (is_right ? size_x - 2 : size_x);
    uint start_j = (is_bottom ? 1 : 0);
    uint end_j = (is_top ? size_y - 2 : size_y);

    for (uint j = start_j; j < end_j; j++)
        for (uint i = start_i; i < end_i; i++)
        {
            vector_cell& center_uv = uv_center.get_cell_ref(i, j);
            vector_cell const center_fg = fg_center.get_cell(i, j);
            scalar_cell const center_p = p_center.get_cell(i, j);
            scalar_cell const right_p = get_neighbor_cell(p_center, p_right, p_right, p_right, p_top, p_top, p_top, p_top, p_top, i, j, RIGHT);
            scalar_cell const top_p = get_neighbor_cell(p_center, p_right, p_right, p_right, p_top, p_top, p_top, p_top, p_top, i, j, TOP);

            center_uv.first = center_fg.first - dt * over_dx * (right_p.value - center_p.value);
            center_uv.second = center_fg.second - dt * over_dy * (top_p.value - center_p.value);
        }

    // compute top strip i = 1, ..., i_max-1, j = j_max for F and set top boundary G = v
    if (is_top)
    {
        uint j = size_y - 2;

        for (uint i = start_i; i < end_i; i++)
        {
            vector_cell& center_uv = uv_center.get_cell_ref(i, j);
            vector_cell const center_fg = fg_center.get_cell(i, j);
            scalar_cell const center_p = p_center.get_cell(i, j);
            scalar_cell const right_p = get_neighbor_cell(p_center, p_right, p_right, p_right, p_top, p_top, p_top, p_top, p_top, i, j, RIGHT);
            scalar_cell const top_p = get_neighbor_cell(p_center, p_right, p_right, p_right, p_top, p_top, p_top, p_top, p_top, i, j, TOP);

            center_uv.first = center_fg.first - dt*over_dx * (right_p.value - center_p.value);
        }
    }

    // compute right strip i = i_max, j = 1, ..., j_max-1 and set right boundary F = u
    if (is_right)
    {
        uint i = size_x - 2;

        for (uint j = start_j; j < end_j; j++)
        {
            vector_cell& center_uv = uv_center.get_cell_ref(i, j);
            vector_cell const center_fg = fg_center.get_cell(i, j);
            scalar_cell const center_p = p_center.get_cell(i, j);
            scalar_cell const right_p = get_neighbor_cell(p_center, p_right, p_right, p_right, p_top, p_top, p_top, p_top, p_top, i, j, RIGHT);
            scalar_cell const top_p = get_neighbor_cell(p_center, p_right, p_right, p_right, p_top, p_top, p_top, p_top, p_top, i, j, TOP);

            center_uv.second = center_fg.second- dt*over_dy * (top_p.value - center_p.value);
        }
    }

    return vector_partition(where, uv_center);
}

vector_partition dispatch_update_velocities(vector_partition const& uv_center, scalar_partition const& center, scalar_partition const& right,
                                        scalar_partition const& top, vector_partition const& fg,
                                        uint global_i, uint global_j, uint i_max, uint j_max, RealType dx, RealType dy, RealType dt)
{
    hpx::shared_future<vector_data> uv_center_data = uv_center.get_data(CENTER);
    hpx::shared_future<scalar_data> center_data = center.get_data(CENTER);
    hpx::shared_future<scalar_data> right_data = right.get_data(RIGHT);
    hpx::shared_future<scalar_data> top_data = top.get_data(TOP);
    hpx::shared_future<vector_data> fg_data = fg.get_data(CENTER);

    hpx::naming::id_type const where = center.get_id();

    return hpx::dataflow(
            hpx::launch::async,
            &update_velocities_action,
            where,
            uv_center_data, center_data, right_data, top_data, fg_data,
            global_i,
            global_j,
            i_max,
            j_max,
            dx,
            dy,
            dt
    );
}

std::pair<RealType, RealType> compute_max_uv_action(hpx::shared_future<vector_data> uv_center_fut)
{
    vector_data center = uv_center_fut.get();

    uint size_x = center.size_x();
    uint size_y = center.size_y();

    std::pair<RealType, RealType> max_uv(0, 0);

    for (uint j = 0; j < size_y; j++)
        for (uint i = 0; i < size_x; i++)
        {
            vector_cell const current = center.get_cell(i, j);
            max_uv.first = (std::abs(current.first) > max_uv.first ? std::abs(current.first) : max_uv.first);
            max_uv.second = (std::abs(current.second) > max_uv.second ? std::abs(current.second) : max_uv.second);
        }

    return max_uv;
}

hpx::future<std::pair<RealType, RealType> > dispatch_compute_max_uv(vector_partition const& uv_center)
{
    hpx::shared_future<vector_data> uv_center_data = uv_center.get_data(CENTER);

    return hpx::dataflow(
            hpx::launch::async,
            &compute_max_uv_action,
            uv_center_data
    );
}

hpx::future<std::pair<RealType, RealType> > custom_grain_size::update_velocities(vector_grid_type& uv_grid,
                                                vector_grid_type const& fg_grid, scalar_grid_type const& p_grid, RealType dt)
{
    for (uint l = 1; l < p.num_partitions_y - 1; l++)
    {
        for (uint k = 1; k < p.num_partitions_x - 1; k++)
        {
            vector_partition& next = uv_grid[get_index(k, l)];

            next =
                hpx::dataflow(
                    hpx::launch::async,
                    &dispatch_update_velocities,
                    next,
                    p_grid[get_index(k, l)], //center
                    p_grid[get_index(k+1, l)], //right
                    p_grid[get_index(k, l+1)], //top
                    fg_grid[get_index(k, l)], //center
                    index[get_index(k, l)].first,
                    index[get_index(k, l)].second,
                    p.i_max,
                    p.j_max,
                    p.dx,
                    p.dy,
                    dt
            );
        }
    }

    hpx::future<std::pair<RealType, RealType> > max_uv = hpx::make_ready_future(std::pair<RealType, RealType>(0, 0));

    for (uint l = 1; l < p.num_partitions_y - 1; l++)
    {
        for (uint k = 1; k < p.num_partitions_x - 1; k++)
        {
            max_uv = hpx::dataflow(
                hpx::launch::async,
                [](hpx::future<std::pair<RealType, RealType> > old_max_uv, hpx::future<std::pair<RealType, RealType> > new_values)
                    -> std::pair<RealType, RealType>
                {
                    std::pair<RealType, RealType> max_uv = old_max_uv.get();
                    std::pair<RealType, RealType> values = new_values.get();

                    max_uv.first = (values.first > max_uv.first ? values.first : max_uv.first);
                    max_uv.second = (values.second > max_uv.second ? values.second : max_uv.second);

                    return max_uv;
                }
                , max_uv
                , hpx::dataflow(
                    hpx::launch::async,
                    &dispatch_compute_max_uv,
                    uv_grid[get_index(k, l)]
                    )
            );
        }
    }

    return max_uv;
}



}//namespace

