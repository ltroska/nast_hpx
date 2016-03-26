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

void with_for_each::set_pressure_on_boundary(scalar_data& p, uint global_i, uint global_j, uint i_max, uint j_max)
{
    uint size_x = p.size_x();
    uint size_y = p.size_y();

    bool is_left = (global_i == 0);
    bool is_right = (global_i + size_x > i_max);
    bool is_bottom = (global_j == 0);
    bool is_top = (global_j + size_y > j_max);

    uint start_i = (is_left ? 1 : 0);
    uint end_i = (is_right ? size_x - 1 : size_x);
    uint start_j = (is_bottom ? 1 : 0);
    uint end_j = (is_top ? size_y - 1 : size_y);

    auto range_j = boost::irange(start_j, end_j);
    auto range_i = boost::irange(start_i, end_i);

    std::vector<hpx::future<void> > futures;

    if (is_left)
    {
        futures.push_back(
            hpx::parallel::for_each(hpx::parallel::par(hpx::parallel::task), boost::begin(range_j), boost::end(range_j),
                [&p](uint j)
                {
                    scalar_cell& cell = p.get_cell_ref(0, j);
                    scalar_cell const cell2 = p.get_cell(1, j);

                    cell.value = cell2.value;
                }
            )
        );
    }

    if (is_right)
    {
        futures.push_back(
            hpx::parallel::for_each(hpx::parallel::par(hpx::parallel::task), boost::begin(range_j), boost::end(range_j),
                [&p, size_x](uint j)
                {
                    scalar_cell& cell = p.get_cell_ref(size_x - 1, j);
                    scalar_cell const cell2 = p.get_cell_ref(size_x - 2, j);

                    cell.value = cell2.value;
                }
            )
        );
    }

    if (is_bottom)
    {
        futures.push_back(
            hpx::parallel::for_each(hpx::parallel::par(hpx::parallel::task), boost::begin(range_i), boost::end(range_i),
                [&p](uint i)
                {
                    scalar_cell& cell = p.get_cell_ref(i, 0);
                    scalar_cell const cell2 = p.get_cell(i, 1);

                    cell.value = cell2.value;
                }
            )
        );
    }

    if (is_top)
    {
       futures.push_back(
            hpx::parallel::for_each(hpx::parallel::par(hpx::parallel::task), boost::begin(range_i), boost::end(range_i),
                [&p, size_y](uint i)
                {
                    scalar_cell& cell = p.get_cell_ref(i, size_y - 1);
                    scalar_cell const cell2 = p.get_cell_ref(i, size_y - 2);

                    cell.value = cell2.value;
                }
            )
        );
    }

    hpx::wait_all(futures);
}

void with_for_each::sor_cycle(scalar_data& p_center, scalar_data const& p_left, scalar_data const& p_right,
                                scalar_data const& p_bottom, scalar_data const& p_top,
                                scalar_data const& rhs_center,
                                uint global_i, uint global_j, uint i_max, uint j_max,
                                RealType omega, RealType dx, RealType dy)
{
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


    RealType dx_sq = std::pow(dx, 2);
    RealType dy_sq = std::pow(dy, 2);
    RealType part1 = 1. - omega;
    RealType part2 = omega * dx_sq * dy_sq / (2. * (dx_sq + dy_sq));

    auto range = boost::irange(0, static_cast<int>(p_center.size()));

  /*  for (uint i = 0; i < size_x; i++)
        for (uint j = 0; j < size_y; j++)
            if (in_range(1, i_max, 1, j_max, global_i + i, global_j + j))
                    {
                        scalar_cell& next_p = p_center.get_cell_ref(i, j);
                        scalar_cell const current_rhs = rhs_center.get_cell(i, j);
                        scalar_cell const left = get_left_neighbor(p_center, p_left, i, j);
                        scalar_cell const right = get_right_neighbor(p_center, p_right, i, j);
                        scalar_cell const bottom = get_bottom_neighbor(p_center, p_bottom, i, j);
                        scalar_cell const top = get_top_neighbor(p_center, p_top, i, j);

                        next_p.value = part1 * next_p.value
                                    + part2 * ( (right.value + left.value) / dx_sq + (top.value + bottom.value) / dy_sq - current_rhs.value);
                    }*/

    hpx::parallel::for_each(hpx::parallel::par, boost::begin(range),
        boost::end(range),
        [&](uint cnt)
        {
            uint const i = cnt%size_x;
            uint const j = cnt/size_x;

            if (in_range(1, i_max, 1, j_max, global_i + i, global_j + j) && (global_i + i + global_j + j)% 2 == 1)
            {
                scalar_cell& next_p = p_center.get_cell_ref(i, j);
                scalar_cell const current_rhs = rhs_center.get_cell(i, j);
                scalar_cell const left = get_left_neighbor(p_center, p_left, i, j);
                scalar_cell const right = get_right_neighbor(p_center, p_right, i, j);
                scalar_cell const bottom = get_bottom_neighbor(p_center, p_bottom, i, j);
                scalar_cell const top = get_top_neighbor(p_center, p_top, i, j);

                next_p.value = part1 * next_p.value
                            + part2 * ( (right.value + left.value) / dx_sq + (top.value + bottom.value) / dy_sq - current_rhs.value);
            }
        }
    );

    hpx::parallel::for_each(hpx::parallel::par, boost::begin(range),
        boost::end(range),
        [&](uint cnt)
        {
            uint const i = cnt%size_x;
            uint const j = cnt/size_x;

            if (in_range(1, i_max, 1, j_max, global_i + i, global_j + j) && (global_i + i + global_j + j)% 2 == 0)
            {
                scalar_cell& next_p = p_center.get_cell_ref(i, j);
                scalar_cell const current_rhs = rhs_center.get_cell(i, j);
                scalar_cell const left = get_left_neighbor(p_center, p_left, i, j);
                scalar_cell const right = get_right_neighbor(p_center, p_right, i, j);
                scalar_cell const bottom = get_bottom_neighbor(p_center, p_bottom, i, j);
                scalar_cell const top = get_top_neighbor(p_center, p_top, i, j);

                next_p.value = part1 * next_p.value
                            + part2 * ( (right.value + left.value) / dx_sq + (top.value + bottom.value) / dy_sq - current_rhs.value);
            }
        }
    );
}


RealType with_for_each::compute_residual(scalar_data const& p_center, scalar_data const& p_left,
                                            scalar_data const& p_right, scalar_data const& p_bottom,
                                            scalar_data const& p_top, scalar_data const& rhs_center,
                                            uint global_i, uint global_j, uint i_max, uint j_max, RealType dx, RealType dy)
{
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

    RealType over_dx_sq = 1./std::pow(dx, 2);
    RealType over_dy_sq = 1./std::pow(dy, 2);

    auto range = boost::irange(0, static_cast<int>(size_x * size_y));

    hpx::future<RealType> local_residual = hpx::parallel::transform_reduce(hpx::parallel::par(hpx::parallel::task), boost::begin(range), boost::end(range),
        [&](uint cnt)
            -> RealType
        {
            uint const i = cnt%size_x;
            uint const j = cnt/size_x;

            if (in_range(start_i, end_i - 1, start_j, end_j - 1, i, j))
            {
                scalar_cell const center = p_center.get_cell(i, j);
                scalar_cell const rhs = rhs_center.get_cell(i, j);
                scalar_cell const left = get_neighbor_cell(p_center, p_left, p_right, p_bottom, p_top, p_top, p_top, p_top, p_top, i, j, LEFT);
                scalar_cell const right = get_neighbor_cell(p_center, p_left, p_right, p_bottom, p_top, p_top, p_top, p_top, p_top, i, j, RIGHT);
                scalar_cell const bottom = get_neighbor_cell(p_center, p_left, p_right, p_bottom, p_top, p_top, p_top, p_top, p_top, i, j, BOTTOM);
                scalar_cell const top = get_neighbor_cell(p_center, p_left, p_right, p_bottom, p_top, p_top, p_top, p_top, p_top, i, j, TOP);

                RealType tmp = (right.value - 2*center.value + left.value)*over_dx_sq + (top.value - 2*center.value + bottom.value)*over_dy_sq - rhs.value;

                return std::pow(tmp, 2);
            }

            return 0.;
        },
        0.,
        [](RealType a, RealType b) -> RealType {return a+b;}
    );

    return local_residual.then(
                                [i_max, j_max](hpx::future<RealType> a) -> RealType
                                    {
                                        return a.get()/(i_max*j_max);
                                    }
                                ).get();
}

void with_for_each::update_velocities(vector_data& uv_center, scalar_data const& p_center, scalar_data const& p_right,
                                        scalar_data const& p_top, vector_data const& fg_center, uint global_i,
                                        uint global_j, uint i_max, uint j_max, RealType dx, RealType dy, RealType dt)
{
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

    auto range = boost::irange(0, static_cast<int>(size_x * size_y));

    hpx::parallel::for_each(hpx::parallel::par, boost::begin(range), boost::end(range),
        [&](uint cnt)
        {
            uint const i = cnt%size_x;
            uint const j = cnt/size_x;

            if (in_range(1, i_max, 1, j_max, global_i + i, global_j + j))
            {
                vector_cell& center_uv = uv_center.get_cell_ref(i, j);
                vector_cell const center_fg = fg_center.get_cell(i, j);
                scalar_cell const center_p = p_center.get_cell(i, j);
                scalar_cell const right_p = get_neighbor_cell(p_center, p_right, p_right, p_right, p_top, p_top, p_top, p_top, p_top, i, j, RIGHT);
                scalar_cell const top_p = get_neighbor_cell(p_center, p_right, p_right, p_right, p_top, p_top, p_top, p_top, p_top, i, j, TOP);

                if (in_range(1, i_max - 1, 1, j_max - 1, global_i + i, global_j + j))
                {
                    center_uv.first = center_fg.first - dt * over_dx * (right_p.value - center_p.value);
                    center_uv.second = center_fg.second - dt * over_dy * (top_p.value - center_p.value);
                }

                if (in_range(1, i_max - 1, j_max, j_max, global_i + i, global_j + j))
                    center_uv.first = center_fg.first - dt*over_dx * (right_p.value - center_p.value);

                if (in_range(i_max, i_max, 1, j_max - 1, global_i + i, global_j + j))
                    center_uv.second = center_fg.second- dt*over_dy * (top_p.value - center_p.value);
            }

        }
    );
}

}//computation

