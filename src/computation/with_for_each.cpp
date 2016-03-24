#include "with_for_each.hpp"

#include <hpx/parallel/algorithm.hpp>
#include <hpx/parallel/algorithms/transform_reduce.hpp>

#include "util/helpers.hpp"
#include "stencils.hpp"

namespace computation {

void with_for_each::set_velocity_on_boundary(vector_data& uv, uint global_i, uint global_j,
                                                uint i_max, uint j_max)
{
    /*
    *@TODO: maybe create new vector_data here
    */
    uint size_x = uv.size_x();
    uint size_y = uv.size_y();

    bool is_left = (global_i == 0);
    bool is_right = (global_i + size_x > i_max);

    bool is_bottom = (global_j == 0);
    bool is_top = (global_j + size_y > j_max);

    uint start_i = (is_left ? 1 : 0);
    uint end_i = (is_right ? size_x - 1 : size_x);
    uint start_j = (is_bottom ? 1 : 0);
    uint end_j = (is_top ? size_y - 1 : size_y);

    auto range_i = boost::irange(start_i, end_i);
    auto range_j = boost::irange(start_j, end_j);

    vector_cell const topright_cell = uv.get_cell(size_x - 2, size_y - 2);

    std::vector<hpx::future<void> > futures;

    if (is_left)
    {
        futures.push_back(
            hpx::parallel::for_each(hpx::parallel::par(hpx::parallel::task), boost::begin(range_j), boost::end(range_j),
                [&uv](uint j)
                {
                    vector_cell& cell = uv.get_cell_ref(0, j);
                    vector_cell const cell2 = uv.get_cell(1, j);

                    cell.first = 0;
                    cell.second = -cell2.second;
                }
            )
        );
    }

    if (is_bottom)
    {
        futures.push_back(
            hpx::parallel::for_each(hpx::parallel::par(hpx::parallel::task), boost::begin(range_i), boost::end(range_i),
                [&uv](uint i)
                {
                    vector_cell& cell = uv.get_cell_ref(i, 0);
                    vector_cell const cell2 = uv.get_cell(i, 1);

                    cell.second = 0;
                    cell.first = -cell2.first;
                }
            )
        );
    }

    if (is_right)
    {
        futures.push_back(
            hpx::parallel::for_each(hpx::parallel::par(hpx::parallel::task), boost::begin(range_j), boost::end(range_j),
                [&uv, size_x](uint j)
                {
                    vector_cell& cell = uv.get_cell_ref(size_x - 2, j);
                    vector_cell& cell2 = uv.get_cell_ref(size_x - 1, j);

                    cell.first = 0;
                    cell2.second = -cell.second;
                }
            )
        );
    }

    if (is_top)
    {
        futures.push_back(
            hpx::parallel::for_each(hpx::parallel::par(hpx::parallel::task), boost::begin(range_i), boost::end(range_i),
                [&uv, size_y](uint i)
                {
                    vector_cell& cell = uv.get_cell_ref(i, size_y - 2);
                    vector_cell& cell2 = uv.get_cell_ref(i, size_y - 1);

                    cell.second = 0;
                    cell2.first = 2. - cell.first;
                }
            )
        );
    }

    hpx::wait_all(futures);

    if (is_right && is_top)
    {
        uv.get_cell_ref(size_x - 2, size_y - 1).first = 2. - topright_cell.first;
        uv.get_cell_ref(size_x - 1, size_y - 2).second = -topright_cell.second;
    }
}

void with_for_each::compute_fg(vector_data& fg, vector_data const& uv_center,
                    vector_data const& uv_left, vector_data const& uv_right,
                    vector_data const& uv_bottom, vector_data const& uv_top,
                    vector_data const& uv_bottomright, vector_data const& uv_topleft,
                    uint global_i, uint global_j, uint i_max, uint j_max, RealType re,
                    RealType dx, RealType dy, RealType dt, RealType alpha)
{

    uint size_x = uv_center.size_x();
    uint size_y = uv_center.size_y();

    bool is_left = (global_i == 0);
    bool is_right = (global_i + size_x > i_max);

    bool is_bottom = (global_j == 0);
    bool is_top = (global_j + size_y > j_max);

    auto range = boost::irange(0, static_cast<int>(size_x * size_y));

    hpx::parallel::for_each(hpx::parallel::par, boost::begin(range),
        boost::end(range),
        [&](uint cnt)
        {
            uint const i = cnt%size_x;
            uint const j = cnt/size_x;

            vector_cell& fg_cell = fg.get_cell_ref(i, j);

            vector_cell const center = uv_center.get_cell(i, j);



            if (in_range(1, i_max - 1, 1, j_max - 1, global_i +i , global_j + j))
            {
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

            if (in_range(1, i_max - 1, j_max, j_max, global_i + i, global_j + j))
            {
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

            if (in_range(i_max, i_max, 1, j_max - 1, global_i + i, global_j + j))
            {
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

            if (in_range(0, 0, 1, j_max, global_i + i, global_j +j))
            {
                fg_cell.first = center.first;
            }

            if (in_range(1, i_max, 0, 0, global_i + i, global_j + j))
            {
                fg_cell.second = center.second;
            }

            if (is_top && i == size_x - 2 && j == size_y - 2)
                fg_cell.second = center.second;

            if (is_right && i == size_x - 2 && j == size_y - 2)
                fg_cell.first = center.first;

        }
    );
}

void with_for_each::compute_rhs(scalar_data& rhs, vector_data const& fg_center, vector_data const& fg_left,
                                    vector_data const& fg_bottom, uint global_i, uint global_j, uint i_max, uint j_max,
                                    RealType dx, RealType dy, RealType dt)
{
    uint size_x = rhs.size_x();
    uint size_y = rhs.size_y();

    auto range = boost::irange(0, static_cast<int>(size_x*size_y));

    hpx::parallel::for_each(hpx::parallel::par, boost::begin(range), boost::end(range),
        [&](uint cnt)
        {
            uint const i = cnt%size_x;
            uint const j = cnt/size_x;

            if (in_range(1, i_max, 1, j_max, global_i + i, global_j + j))
            {
                vector_cell const center = fg_center.get_cell(i, j);
                vector_cell const left = get_neighbor_cell(fg_center, fg_left, fg_left, fg_bottom, fg_bottom, fg_bottom, fg_bottom, fg_bottom, fg_bottom, i, j, LEFT);
                vector_cell const bottom = get_neighbor_cell(fg_center, fg_left, fg_left, fg_bottom, fg_bottom, fg_bottom, fg_bottom, fg_bottom, fg_bottom, i, j, BOTTOM);

                rhs.get_cell_ref(i, j).value = 1./dt * ( (center.first - left.first)/dx + (center.second - bottom.second)/dy);
            }
        }
    );
}

}//computation

