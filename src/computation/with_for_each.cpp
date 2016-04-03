#include "with_for_each.hpp"

#include <hpx/parallel/algorithm.hpp>
#include <hpx/parallel/algorithms/transform_reduce.hpp>
#include <hpx/parallel/executors/static_chunk_size.hpp>

#include "util/helpers.hpp"
#include "stencils.hpp"

namespace computation {

void with_for_each::set_boundary(vector_data& uv_center, vector_data const& uv_left, vector_data const& uv_right, vector_data const& uv_bottom, vector_data const& uv_top,
                                    scalar_data& temperature, std::vector<std::bitset<5> > const& flag_data, boundary_data const& data_type,
                                    boundary_data const& temp_data_type, boundary_data const& u_bnd, boundary_data const& v_bnd, boundary_data const& temp_bnd,
                                    uint global_i, uint global_j, uint i_max, uint j_max, RealType dx, RealType dy)
{
    /*
    *@TODO: maybe create new vector_data here
    */
    uint size_x = uv_center.size_x();
    uint size_y = uv_center.size_y();

    bool is_left = (global_i == 0);
    bool is_right = (global_i + size_x > i_max);

    bool is_bottom = (global_j == 0);
    bool is_top = (global_j + size_y > j_max);

    uint start_i = (is_left ? 1 : 0);
    uint end_i = (is_right ? size_x - 1 : size_x);
    uint start_j = (is_bottom ? 1 : 0);
    uint end_j = (is_top ? size_y - 1 : size_y);

    auto range = boost::irange(0, static_cast<int>(size_x*size_y));

    vector_cell const topright_cell = uv_center.get_cell(size_x - 2, size_y - 2);

    hpx::parallel::for_each(hpx::parallel::par, boost::begin(range), boost::end(range),
        [&](uint cnt)
        {
            uint const i = cnt%size_x;
            uint const j = cnt/size_x;

            std::bitset<5> cell_type = flag_data[j*size_x + i];

            //obstacle
            if (in_range(1, i_max, 1, j_max, global_i + i, global_j + j))
            {
                //north
                if (cell_type == std::bitset<5>("00001"))
                {
                    vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                    vector_cell const top_cell = get_top_neighbor(uv_center, uv_top, i, j);

                    curr_cell.second = 0;
                    curr_cell.first = -top_cell.first;
                }

                //east
                if (cell_type == std::bitset<5>("01000"))
                {
                    vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                    vector_cell const right_cell = get_right_neighbor(uv_center, uv_right, i, j);

                    curr_cell.first = 0;
                    curr_cell.second = -right_cell.second;
                }

                //south
                if (cell_type == std::bitset<5>("00010"))
                {
                    vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                    vector_cell const bottom_cell = get_bottom_neighbor(uv_center, uv_bottom, i, j);

                    curr_cell.first = -bottom_cell.first;
                }

                if (cell_type == std::bitset<5>("11110") && global_j + j != j_max)
                {
                    vector_cell& curr_cell = uv_center.get_cell_ref(i, j);

                    curr_cell.second = 0;
                }

                //west
                if (cell_type == std::bitset<5>("00100"))
                {
                    vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                    vector_cell const left_cell = get_left_neighbor(uv_center, uv_left, i, j);

                    curr_cell.second = -left_cell.second;
                }

                if (cell_type == std::bitset<5>("10111") && global_i + i != i_max)
                {
                    vector_cell& curr_cell = uv_center.get_cell_ref(i, j);

                    curr_cell.first = 0;
                }

                //north east
                if (cell_type == std::bitset<5>("01001"))
                {
                    vector_cell& curr_cell = uv_center.get_cell_ref(i, j);

                    curr_cell.first = 0;
                    curr_cell.second = 0;
                }

                //south east
                if (cell_type == std::bitset<5>("01010"))
                {
                    vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                    vector_cell const right_cell = get_right_neighbor(uv_center, uv_right, i, j);

                    curr_cell.first = 0;
                    curr_cell.second = -right_cell.second;
                }

                //south west
                if (cell_type == std::bitset<5>("00110"))
                {
                    vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                    vector_cell const bottom_cell = get_bottom_neighbor(uv_center, uv_bottom, i, j);
                    vector_cell const left_cell = get_left_neighbor(uv_center, uv_left, i, j);

                    curr_cell.first = -bottom_cell.first;
                    curr_cell.second = -left_cell.second;
                }

                //north west
                if (cell_type == std::bitset<5>("00101"))
                {
                    vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                    vector_cell const top_cell = get_top_neighbor(uv_center, uv_top, i, j);

                    curr_cell.first = -top_cell.first;
                    curr_cell.second = 0;
                }
            }
            //boundary
            else
            {
                //left
                if (in_range(0, 0, 1, j_max, global_i + i, global_j + j))
                {
                    if (data_type.left == 1)
                    {
                        vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                        vector_cell const right_cell = get_right_neighbor(uv_center, uv_right, i, j);

                        curr_cell.first = 0;
                        curr_cell.second = -right_cell.second;
                    }

                    if (data_type.left == 2)
                    {
                        vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                        vector_cell const right_cell = get_right_neighbor(uv_center, uv_right, i, j);

                        curr_cell.first = 0;
                        curr_cell.second = right_cell.second;
                    }

                    if (data_type.left == 3)
                    {
                        vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                        vector_cell const right_cell = get_right_neighbor(uv_center, uv_right, i, j);

                        curr_cell.first = right_cell.first;
                        curr_cell.second = right_cell.second;
                    }

                    if (data_type.left == 4)
                    {
                        vector_cell& curr_cell = uv_center.get_cell_ref(i, j);

                        curr_cell.first = u_bnd.left;
                        curr_cell.second = v_bnd.left;
                    }

                    if (temp_data_type.left != -1)
                    {
                        //Dirichlet
                        if (temp_data_type.left == 1)
                        {
                            scalar_cell& curr_cell = temperature.get_cell_ref(i, j);
                            scalar_cell const right_cell = temperature.get_cell(i+1, j);

                            curr_cell.value = 2*temp_bnd.left*((j - 0.5)*dy) - right_cell.value;
                        }

                        //Neumann
                        else if (temp_data_type.left == 2)
                        {
                            scalar_cell& curr_cell = temperature.get_cell_ref(i, j);
                            scalar_cell const right_cell = temperature.get_cell(i+1, j);

                            curr_cell.value = right_cell.value + dx*temp_bnd.left*((j-0.5)*dy);
                        }
                    }
                }
                //right
                if (in_range(i_max + 1, i_max + 1, 1, j_max, global_i + i, global_j + j))
                {
                    if (data_type.right == 1)
                    {
                        vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                        vector_cell& left_cell = uv_center.get_cell_ref(i - 1, j);

                        left_cell.first = 0;
                        curr_cell.second = -left_cell.second;
                    }

                    if (data_type.right == 2)
                    {
                        vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                        vector_cell& left_cell = uv_center.get_cell_ref(i - 1, j);

                        left_cell.first = 0;
                        curr_cell.second = left_cell.second;
                    }

                    if (data_type.right == 3)
                    {
                        vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                        vector_cell& left_cell = uv_center.get_cell_ref(i - 1, j);
                        vector_cell const left2_cell = get_left_neighbor(uv_center, uv_left, i-1, j);

                        left_cell.first = left2_cell.first;
                        curr_cell.second = left_cell.second;
                    }

                    if (data_type.right == 4)
                    {
                        vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                        vector_cell& left_cell = uv_center.get_cell_ref(i - 1, j);

                        left_cell.first = u_bnd.right;
                        curr_cell.second = v_bnd.right;
                    }

                    if (temp_data_type.right != -1)
                    {
                        //Dirichlet
                        if (temp_data_type.right == 1)
                        {
                            scalar_cell& curr_cell = temperature.get_cell_ref(i, j);
                            scalar_cell const left_cell = temperature.get_cell(i-1, j);

                            curr_cell.value = 2*temp_bnd.right * ((j - 0.5)*dy) - left_cell.value;
                        }

                        //Neumann
                        else if (temp_data_type.right == 2)
                        {
                            scalar_cell& curr_cell = temperature.get_cell_ref(i, j);
                            scalar_cell const left_cell = temperature.get_cell(i-1, j);

                            curr_cell.value = left_cell.value + dx*temp_bnd.right*((j-0.5)*dy);
                        }
                    }
                }

                //bottom
                if (in_range(1, i_max, 0, 0, global_i + i, global_j + j))
                {
                    if (data_type.bottom == 1)
                    {
                        vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                        vector_cell const top_cell = get_top_neighbor(uv_center, uv_top, i, j);

                        curr_cell.first = -top_cell.first;
                        curr_cell.second = 0;
                    }

                    if (data_type.bottom == 2)
                    {
                        vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                        vector_cell const top_cell = get_top_neighbor(uv_center, uv_top, i, j);

                        curr_cell.first = top_cell.first;
                        curr_cell.second = 0;
                    }

                    if (data_type.bottom == 3)
                    {
                        vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                        vector_cell const top_cell = get_top_neighbor(uv_center, uv_top, i, j);

                        curr_cell.first = top_cell.first;
                        curr_cell.second = top_cell.second;
                    }

                    if (data_type.bottom == 4)
                    {
                        vector_cell& curr_cell = uv_center.get_cell_ref(i, j);

                        curr_cell.first = u_bnd.bottom;
                        curr_cell.second = v_bnd.bottom;
                    }

                    if (temp_data_type.bottom != -1)
                    {
                        //Dirichlet
                        if (temp_data_type.bottom == 1)
                        {
                            scalar_cell& curr_cell = temperature.get_cell_ref(i, j);
                            scalar_cell const top_cell = temperature.get_cell(i, j+1);

                            curr_cell.value = 2*temp_bnd.bottom * ((i - 0.5)*dx) - top_cell.value;
                        }

                        //Neumann
                        else if (temp_data_type.bottom == 2)
                        {
                            scalar_cell& curr_cell = temperature.get_cell_ref(i, j);
                            scalar_cell const top_cell = temperature.get_cell(i, j+1);

                            curr_cell.value = top_cell.value + dy*temp_bnd.bottom*((i-0.5)*dx);
                        }
                    }
                }

                //top
                if (in_range(1, i_max, j_max + 1, j_max + 1, global_i + i, global_j + j))
                {
                    if (data_type.top == 1)
                    {
                        vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                        vector_cell& bottom_cell = uv_center.get_cell_ref(i, j - 1);

                        curr_cell.first = 2*u_bnd.top - bottom_cell.first;
                        bottom_cell.second = 0;
                    }

                    if (data_type.top == 2)
                    {
                        vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                        vector_cell& bottom_cell = uv_center.get_cell_ref(i, j - 1);

                        curr_cell.first = bottom_cell.first;
                        bottom_cell.second = 0;
                    }

                    if (data_type.top == 3)
                    {
                        vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                        vector_cell& bottom_cell = uv_center.get_cell_ref(i, j - 1);
                        vector_cell const bottom2_cell = get_bottom_neighbor(uv_center, uv_bottom, i, j-1);

                        curr_cell.first = bottom_cell.first;
                        bottom_cell.second = bottom2_cell.second;
                    }

                    if (data_type.top == 4)
                    {
                        vector_cell& curr_cell = uv_center.get_cell_ref(i, j);
                        vector_cell& bottom_cell = uv_center.get_cell_ref(i, j - 1);

                        curr_cell.first = u_bnd.top;
                        bottom_cell.second = v_bnd.top;
                    }

                    if (temp_data_type.top != -1)
                    {
                        //Dirichlet
                        if (temp_data_type.top == 1)
                        {
                            scalar_cell& curr_cell = temperature.get_cell_ref(i, j);
                            scalar_cell const bottom_cell = temperature.get_cell(i, j-1);

                            curr_cell.value = 2*temp_bnd.top * ((i - 0.5)*dx) - bottom_cell.value;
                        }

                        //Neumann
                        else if (temp_data_type.top == 2)
                        {
                            scalar_cell& curr_cell = temperature.get_cell_ref(i, j);
                            scalar_cell const bottom_cell = temperature.get_cell(i, j-1);

                            curr_cell.value = bottom_cell.value + dy*temp_bnd.top*((i-0.5)*dx);
                        }
                    }
                }
            }
        }
    );

    if (is_right && is_top)
    {
        if (data_type.top == 1)
            uv_center.get_cell_ref(size_x - 2, size_y - 1).first = 2*u_bnd.top - topright_cell.first;

        if (data_type.top == 2 || data_type.top == 3)
            uv_center.get_cell_ref(size_x - 2, size_y - 1).first = topright_cell.first;

        if (data_type.right == 1)
            uv_center.get_cell_ref(size_x - 1, size_y - 2).second = -topright_cell.second;

        if (data_type.right == 2 || data_type.right == 3)
            uv_center.get_cell_ref(size_x - 1, size_y - 2).second = topright_cell.second;
    }

}

void with_for_each::compute_fg(vector_data& fg, vector_data const& uv_center,
                        vector_data const& uv_left, vector_data const& uv_right,
                        vector_data const& uv_bottom, vector_data const& uv_top,
                        vector_data const& uv_bottomright, vector_data const& uv_topleft,
                        scalar_data const& temp_center, scalar_data const& temp_right,
                        scalar_data const& temp_top,
                        std::vector<std::bitset<5> > const& flag_data,
                        uint global_i, uint global_j, uint i_max, uint j_max, RealType re,
                        RealType gx, RealType gy, RealType beta,
                        RealType dx, RealType dy, RealType dt, RealType alpha)
{

    uint size_x = uv_center.size_x();
    uint size_y = uv_center.size_y();

    bool is_left = (global_i == 0);
    bool is_right = (global_i + size_x > i_max);

    bool is_bottom = (global_j == 0);
    bool is_top = (global_j + size_y > j_max);

    auto range = boost::irange(0, static_cast<int>(size_x * size_y));

    hpx::parallel::for_each(hpx::parallel::par, boost::begin(range), boost::end(range),
        [&](uint cnt)
        {
            uint const i = cnt%size_x;
            uint const j = cnt/size_x;

            vector_cell& fg_cell = fg.get_cell_ref(i, j);

            vector_cell const center = uv_center.get_cell(i, j);

            std::bitset<5> cell_type = flag_data[j*size_x + i];

            //north
            if (!cell_type.test(4) && cell_type.test(0))
            {
                vector_cell& curr_cell = fg.get_cell_ref(i, j);
                vector_cell const curr_uv = uv_center.get_cell(i, j);

                curr_cell.second = curr_uv.second;
            }

            //east
            if (!cell_type.test(4) && cell_type.test(3))
            {
                vector_cell& curr_cell = fg.get_cell_ref(i, j);
                vector_cell const curr_uv = uv_center.get_cell(i, j);

                curr_cell.first = curr_uv.first;
            }

            if (cell_type.test(4))
            {
                vector_cell const left = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, LEFT);
                vector_cell const right = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, RIGHT);
                vector_cell const bottom = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, BOTTOM);
                vector_cell const top = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, TOP);
                vector_cell const bottomright = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, BOTTOM_RIGHT);
                vector_cell const topleft = get_neighbor_cell(uv_center, uv_left, uv_right, uv_bottom, uv_top, uv_top, uv_bottomright, uv_topleft, uv_topleft, i, j, TOP_LEFT);

                scalar_cell const t_center = temp_center.get_cell(i, j);
                scalar_cell const t_right = get_right_neighbor(temp_center, temp_right, i, j);
                scalar_cell const t_top = get_top_neighbor(temp_center, temp_top, i, j);

                if (!cell_type.test(3))
                {
                    fg_cell.first = center.first;
                }
                else
                {
                    fg_cell.first = center.first + dt * (
                                    1./re * (second_derivative_fwd_bkwd_x(right.first, center.first, left.first, dx)
                                                + second_derivative_fwd_bkwd_y(top.first, center.first, bottom.first, dy))

                                    - first_derivative_of_square_x(right.first, center.first, left.first, dx, alpha)
                                    - first_derivative_of_product_y(right.second, center.second, bottom.second, bottomright.second,
                                                                    bottom.first, center.first, top.first, dy, alpha)

                                    + gx
                                    )
                                    - beta * dt / 2. * (t_center.value + t_right.value) * gx;
                }

                if (!cell_type.test(0))
                {
                    fg_cell.second = center.second;
                }
                else
                {
                    fg_cell.second = center.second + dt * (
                                    1./re * (second_derivative_fwd_bkwd_x(right.second, center.second, left.second, dx)
                                                + second_derivative_fwd_bkwd_y(top.second, center.second, bottom.second, dy))

                                    - first_derivative_of_product_x(left.first, center.first, top.first, topleft.first,
                                                                    left.second, center.second, right.second, dx, alpha)
                                    - first_derivative_of_square_y(top.second, center.second, bottom.second, dy, alpha)

                                    + gy
                                    )
                                    - beta * dt / 2. * (t_center.value + t_top.value) * gy;
                }
            }

        }
    );
}

void with_for_each::compute_temp(scalar_data& temp_center, scalar_data const& temp_left, scalar_data const& temp_right,
                                    scalar_data const& temp_bottom, scalar_data const& temp_top,
                                    vector_data const& uv_center, vector_data const& uv_left,
                                    vector_data const& uv_bottom,
                                    uint global_i, uint global_j, uint i_max, uint j_max, RealType re, RealType pr,
                                    RealType dx, RealType dy, RealType dt, RealType alpha)
{
    uint size_x = uv_center.size_x();
    uint size_y = uv_center.size_y();

    bool is_left = (global_i == 0);
    bool is_right = (global_i + size_x > i_max);

    bool is_bottom = (global_j == 0);
    bool is_top = (global_j + size_y > j_max);

    auto range = boost::irange(0, static_cast<int>(size_x * size_y));

    hpx::parallel::for_each(hpx::parallel::par, boost::begin(range), boost::end(range),
        [&](uint cnt)
        {
            uint const i = cnt%size_x;
            uint const j = cnt/size_x;

            scalar_cell& temp_cell = temp_center.get_cell_ref(i, j);
            scalar_cell const temp_left_cell = get_left_neighbor(temp_center, temp_left, i, j);
            scalar_cell const temp_right_cell = get_right_neighbor(temp_center, temp_right, i, j);
            scalar_cell const temp_bottom_cell = get_bottom_neighbor(temp_center, temp_bottom, i, j);
            scalar_cell const temp_top_cell = get_top_neighbor(temp_center, temp_top, i, j);

            vector_cell const uv_center_cell = uv_center.get_cell(i, j);
            vector_cell const uv_left_cell = get_left_neighbor(uv_center, uv_left, i, j);
            vector_cell const uv_bottom_cell = get_bottom_neighbor(uv_center, uv_bottom, i, j);

            if (in_range(1, i_max, 1, j_max, global_i + i, global_j + j))
            {
                temp_cell.value = dt * (1./re*1./pr * (second_derivative_fwd_bkwd_x(temp_right_cell.value, temp_cell.value, temp_left_cell.value, dx)
                                                        + second_derivative_fwd_bkwd_y(temp_top_cell.value, temp_cell.value, temp_bottom_cell.value, dy)
                                                        )
                                       // + (i == 1 ? 1 : 0)
                                        - first_derivative_u_temp_x(uv_center_cell.first, uv_left_cell.first, temp_cell.value, temp_right_cell.value, temp_left_cell.value, dx, alpha)
                                        - first_derivative_v_temp_y(uv_center_cell.second, uv_bottom_cell.second, temp_cell.value, temp_top_cell.value, temp_bottom_cell.value, dy, alpha)
                                        )
                                        + temp_cell.value;
            }
        }
    );
}

void with_for_each::compute_rhs(scalar_data& rhs, vector_data const& fg_center, vector_data const& fg_left,
                                    vector_data const& fg_bottom, std::vector<std::bitset<5> > const& flag_data, uint global_i, uint global_j, uint i_max, uint j_max,
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

            std::bitset<5> cell_type = flag_data[j*size_x + i];

            if (cell_type.test(4))
            {
                vector_cell const center = fg_center.get_cell(i, j);
                vector_cell const left = get_neighbor_cell(fg_center, fg_left, fg_left, fg_bottom, fg_bottom, fg_bottom, fg_bottom, fg_bottom, fg_bottom, i, j, LEFT);
                vector_cell const bottom = get_neighbor_cell(fg_center, fg_left, fg_left, fg_bottom, fg_bottom, fg_bottom, fg_bottom, fg_bottom, fg_bottom, i, j, BOTTOM);

                rhs.get_cell_ref(i, j).value = 1./dt * ( (center.first - left.first)/dx + (center.second - bottom.second)/dy);
            }
        }
    );
}

void with_for_each::set_pressure_on_boundary(scalar_data& p_center, scalar_data const& p_left, scalar_data const& p_right,
                                                scalar_data const& p_bottom, scalar_data const& p_top, std::vector<std::bitset<5> > const& flag_data,
                                                uint global_i, uint global_j, uint i_max, uint j_max)
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

    auto range = boost::irange(0, static_cast<int>(size_x*size_y));

    hpx::parallel::for_each(hpx::parallel::par, boost::begin(range), boost::end(range),
        [&](uint cnt)
        {
            uint const i = cnt%size_x;
            uint const j = cnt/size_x;

            scalar_cell& curr_cell = p_center.get_cell_ref(i, j);

            std::bitset<5> cell_type = flag_data[j*size_x + i];

            //east
            if (cell_type == std::bitset<5>("01000"))
                curr_cell.value = get_right_neighbor(p_center, p_right, i, j).value;

            //west
            else if (cell_type == std::bitset<5>("00100"))
                curr_cell.value = get_left_neighbor(p_center, p_left, i, j).value;

            //south
            else if (cell_type == std::bitset<5>("00010"))
                curr_cell.value = get_bottom_neighbor(p_center, p_bottom, i, j).value;

            //north
            else if (cell_type == std::bitset<5>("00001"))
                curr_cell.value = get_top_neighbor(p_center, p_top, i, j).value;

            //NE
            else if (cell_type == std::bitset<5>("01001"))
                curr_cell.value = (get_top_neighbor(p_center, p_top, i, j).value + get_right_neighbor(p_center, p_right, i, j).value) / 2.;

            //SE
            else if (cell_type == std::bitset<5>("01010"))
                curr_cell.value = (get_bottom_neighbor(p_center, p_bottom, i, j).value + get_right_neighbor(p_center, p_right, i, j).value) / 2.;

            //SW
            else if (cell_type == std::bitset<5>("00110"))
                curr_cell.value = (get_bottom_neighbor(p_center, p_bottom, i, j).value + get_left_neighbor(p_center, p_left, i, j).value) / 2.;

            //NW
            else if (cell_type == std::bitset<5>("00101"))
                curr_cell.value = (get_top_neighbor(p_center, p_top, i, j).value + get_left_neighbor(p_center, p_left, i, j).value) / 2.;

        }
    );
}

void with_for_each::sor_cycle(scalar_data& p_center, scalar_data const& p_left, scalar_data const& p_right,
                                scalar_data const& p_bottom, scalar_data const& p_top,
                                scalar_data const& rhs_center, std::vector<std::bitset<5> > const& flag_data,
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

    for (uint i = 0; i < size_x; i++)
        for (uint j = 0; j < size_y; j++)
        {
            std::bitset<5> cell_type = flag_data[j*size_x + i];

            if (cell_type.test(4))
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

  /*  hpx::parallel::for_each(hpx::parallel::par, boost::begin(range),
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
    );*/
}


RealType with_for_each::compute_residual(scalar_data const& p_center, scalar_data const& p_left,
                                            scalar_data const& p_right, scalar_data const& p_bottom,
                                            scalar_data const& p_top, scalar_data const& rhs_center,
                                            std::vector<std::bitset<5> > const& flag_data,
                                            uint global_i, uint global_j, uint i_max, uint j_max, RealType dx,
                                            RealType dy)
{
    uint size_x = p_center.size_x();
    uint size_y = p_center.size_y();

    bool is_left = (global_i == 0);
    bool is_right = (global_i + size_x > i_max);
    bool is_bottom = (global_j == 0);
    bool is_top = (global_j + size_y > j_max);

    RealType over_dx_sq = 1./std::pow(dx, 2);
    RealType over_dy_sq = 1./std::pow(dy, 2);

    auto range = boost::irange(0, static_cast<int>(size_x * size_y));

    hpx::future<RealType> local_residual = hpx::parallel::transform_reduce(hpx::parallel::par(hpx::parallel::task), boost::begin(range), boost::end(range),
        [&](uint cnt)
            -> RealType
        {
            uint const i = cnt%size_x;
            uint const j = cnt/size_x;

            std::bitset<5> cell_type = flag_data[j*size_x + i];

            if (cell_type.test(4))
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
                                        scalar_data const& p_top, vector_data const& fg_center, std::vector<std::bitset<5> > const& flag_data,
                                        uint global_i, uint global_j, uint i_max, uint j_max, RealType dx, RealType dy, RealType dt)
{
    uint size_x = uv_center.size_x();
    uint size_y = uv_center.size_y();

    bool is_left = (global_i == 0);
    bool is_right = (global_i + size_x > i_max);
    bool is_bottom = (global_j == 0);
    bool is_top = (global_j + size_y > j_max);

    RealType over_dx = 1./dx;
    RealType over_dy = 1./dy;

    auto range = boost::irange(0, static_cast<int>(size_x * size_y));

    hpx::parallel::for_each(hpx::parallel::par, boost::begin(range), boost::end(range),
        [&](uint cnt)
        {
            uint const i = cnt%size_x;
            uint const j = cnt/size_x;

            std::bitset<5> cell_type = flag_data[j*size_x + i];

            if (cell_type.test(4))
            {
                vector_cell& center_uv = uv_center.get_cell_ref(i, j);
                vector_cell const center_fg = fg_center.get_cell(i, j);
                scalar_cell const center_p = p_center.get_cell(i, j);
                scalar_cell const right_p = get_neighbor_cell(p_center, p_right, p_right, p_right, p_top, p_top, p_top, p_top, p_top, i, j, RIGHT);
                scalar_cell const top_p = get_neighbor_cell(p_center, p_right, p_right, p_right, p_top, p_top, p_top, p_top, p_top, i, j, TOP);

                if (cell_type.test(3))
                    center_uv.first = center_fg.first - dt * over_dx * (right_p.value - center_p.value);

                if (cell_type.test(0))
                    center_uv.second = center_fg.second - dt * over_dy * (top_p.value - center_p.value);
            }

        }
    );
}

void with_for_each::compute_stream_vorticity_heat(scalar_data& stream_center, scalar_data& vorticity_center, scalar_data& heat_center,
                                                    scalar_data const& stream_bottom, scalar_data const& heat_bottom,
                                                    vector_data const& uv_center, vector_data const& uv_right, vector_data const& uv_top,
                                                    scalar_data const& temp_center, scalar_data const& temp_right,
                                                    std::vector<std::bitset<5> > const& flag_data,
                                                    uint global_i, uint global_j, uint i_max, uint j_max, RealType re, RealType pr,
                                                    RealType dx, RealType dy)
{
    uint size_x = stream_center.size_x();
    uint size_y = stream_center.size_y();

    bool is_left = (global_i == 0);
    bool is_right = (global_i + size_x > i_max);
    bool is_bottom = (global_j == 0);
    bool is_top = (global_j + size_y > j_max);

    uint start_i = (is_left ? 1 : 0);
    uint end_i = (is_right ? size_x - 1 : size_x);
    uint start_j = (is_bottom ? 1 : 0);
    uint end_j = (is_top ? size_y - 1 : size_y);

    auto range = boost::irange(0, static_cast<int>(size_x*size_y));

    for (uint j = 0; j < size_y; j++)
        for (uint i = 0; i < size_x; i++)
        {
            std::bitset<5> cell_type = flag_data[j*size_x + i];


            if (in_range(0, i_max, 1, j_max, global_i + i, global_j + j))
            {
                if (cell_type.test(4))
                {
                    stream_center.get_cell_ref(i, j).value = get_bottom_neighbor(stream_center, stream_bottom, i, j).value + uv_center.get_cell(i, j).first*dy;
                    heat_center.get_cell_ref(i, j).value = get_bottom_neighbor(heat_center, heat_bottom, i, j).value
                                                            + dy * (re * pr * uv_center.get_cell(i, j).first
                                                                        * (get_right_neighbor(temp_center, temp_right, i, j).value + temp_center.get_cell(i, j).value) / 2.
                                                                        - (get_right_neighbor(temp_center, temp_right, i, j).value - temp_center.get_cell(i, j).value) / dx
                                                                    );
                }
                else
                {
                     stream_center.get_cell_ref(i, j).value = get_bottom_neighbor(stream_center, stream_bottom, i, j).value;
                     heat_center.get_cell_ref(i, j).value = get_bottom_neighbor(heat_center, heat_bottom, i, j).value;
                }
            }

            if (in_range(1, i_max - 1, 1, j_max - 1, global_i + i, global_j + j))
            {
                if (cell_type.test(4))
                {
                    vector_cell const curr_cell = uv_center.get_cell(i, j);
                    vorticity_center.get_cell_ref(i, j).value = (get_top_neighbor(uv_center, uv_top, i, j).first - curr_cell.first)/dy
                                                                - (get_right_neighbor(uv_center, uv_right, i, j).second - curr_cell.second)/dx;
                }
                else
                {
                    vorticity_center.get_cell_ref(i, j).value = 0;
                }
            }
        }


}

}//computation

