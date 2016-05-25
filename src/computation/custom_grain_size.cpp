#include "custom_grain_size.hpp"

#include "util/helpers.hpp"
#include "cell_operations.hpp"

//TODO: only use get_*_neighbor where necessary

namespace computation {

vector_partition custom_grain_size::set_velocity_for_boundary_and_obstacles(
    vector_partition const& middle, vector_partition const& left,
    vector_partition const& right, vector_partition const& bottom,
    vector_partition const& top,
    std::vector<std::vector<std::pair<uint, uint> > > const& boundary,
    std::vector<std::pair<uint, uint> > const& obstacle,
    std::vector<std::pair<uint, uint> > const& fluid,
    std::vector<std::bitset<5> > const& flag_data,
    boundary_data const& type, boundary_data const& u, boundary_data const& v)
{
    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [middle, flag_data, boundary, obstacle, fluid, type, u, v]
            (vector_data m, vector_data const& l, vector_data const& r,
                vector_data const& b, vector_data const& t)
            -> vector_partition
            {
                auto size_x = m.size_x();
                auto size_y = m.size_y();

                for (auto& idx_pair : obstacle)
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    set_velocity_for_obstacle(
                        m(i, j),
                        get_left_neighbor(m, l, i, j),
                        get_right_neighbor(m, r, i, j),
                        get_bottom_neighbor(m, b, i, j),
                        get_top_neighbor(m, t, i, j),
                        flag_data[j * size_x + i]);
                }

                for (auto& idx_pair : fluid)
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    set_velocity_for_fluid(
                        m(i, j),
                        get_left_neighbor(m, l, i, j),
                        get_right_neighbor(m, r, i, j),
                        get_bottom_neighbor(m, b, i, j),
                        get_top_neighbor(m, t, i, j),
                        flag_data[j * size_x + i]);
                }


                //left
                for (auto& idx_pair : boundary[0])
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    set_velocity_for_left_boundary(
                        m(i, j),
                        m(i + 1, j),
                        type.left, u.left, v.left);
                }

                //right
                for (auto& idx_pair : boundary[1])
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    set_velocity_for_right_boundary(
                        m(i, j),
                        m(i - 1, j),
                        get_left_neighbor(m, l, i - 1, j),
                        type.right, u.right, v.right);
                }

                //bottom
                for (auto& idx_pair : boundary[2])
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    set_velocity_for_bottom_boundary(
                        m(i, j),
                        m(i, j + 1),
                        type.bottom, u.bottom, v.bottom);
                }

                //top
                for (auto& idx_pair : boundary[3])
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    set_velocity_for_top_boundary(
                        m(i, j),
                        m(i, j - 1),
                        get_top_neighbor(m, r, i, j - 1),
                        type.top, u.top, v.top);
                }

                return middle;
            }
        ),
        middle.get_data(CENTER),
        left.get_data(LEFT),
        right.get_data(RIGHT),
        bottom.get_data(BOTTOM),
        top.get_data(TOP)
    );
}

scalar_partition custom_grain_size::set_temperature_for_boundary_and_obstacles(
    scalar_partition const& middle, scalar_partition const& left,
    scalar_partition const& right, scalar_partition const& bottom,
    scalar_partition const& top,
    std::vector<std::vector<std::pair<uint, uint> > > const& boundary,
    boundary_data const& boundary_data_type,
    boundary_data const& temperature_boundary_data,
    uint global_i, uint global_j, RealType dx, RealType dy)
{
    //TODO: add temperature for obstacle cells
    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [middle, boundary_data_type, temperature_boundary_data,
             boundary, global_i, global_j, dx, dy]
            (scalar_data next, scalar_data const& l,
                scalar_data const& r, scalar_data const& b,
                scalar_data const& t)
            -> scalar_partition
            {
                uint size_x = next.size_x();
                uint size_y = next.size_y();

                //left
                for (auto& idx_pair : boundary[0])
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    set_temperature_for_boundary(
                        next(i, j),
                        next(i + 1, j),
                        boundary_data_type.left,
                        temperature_boundary_data.left, j, dx, dy);
                }

                //right
                for (auto& idx_pair : boundary[1])
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    set_temperature_for_boundary(
                        next(i, j),
                        next(i - 1, j),
                        boundary_data_type.right,
                        temperature_boundary_data.right, j, dx, dy);
                }

                //bottom
                for (auto& idx_pair : boundary[2])
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    set_temperature_for_boundary(
                        next(i, j),
                        next(i, j + 1),
                        boundary_data_type.bottom,
                        temperature_boundary_data.bottom, i, dy, dx);
                }

                //top
                for (auto& idx_pair : boundary[3])
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    set_temperature_for_boundary(
                        next(i, j),
                        next(i, j - 1),
                        boundary_data_type.top,
                        temperature_boundary_data.top, i, dy, dx);
                }

                return middle;
            }
        ),
        middle.get_data(CENTER),
        left.get_data(LEFT),
        right.get_data(RIGHT),
        bottom.get_data(BOTTOM),
        top.get_data(TOP)
    );
}

vector_partition custom_grain_size::compute_fg_on_fluid_cells(
    vector_partition const& middle_fg,
    vector_partition const& middle_uv, vector_partition const& left_uv,
    vector_partition const& right_uv, vector_partition const& bottom_uv,
    vector_partition const& top_uv, vector_partition const& bottomright_uv,
    vector_partition const& topleft_uv,
    scalar_partition const& middle_temperature,
    scalar_partition const& right_temperature,
    scalar_partition const& top_temperature,
    std::vector<std::vector<std::pair<uint, uint> > > const& boundary,
    std::vector<std::pair<uint, uint> > const& obstacle,
    std::vector<std::pair<uint, uint> > const& fluid,
    std::vector<std::bitset<5> > const& flag_data, RealType re, RealType gx,
    RealType gy, RealType beta, RealType dx, RealType dy, RealType dt,
    RealType alpha)
{
    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [middle_fg, flag_data, boundary, obstacle, fluid, re, gx, gy, beta,
                dx, dy, dt, alpha]
            (vector_data m_fg, vector_data const& m_uv, vector_data const& l_uv,
                vector_data const& r_uv, vector_data const& b_uv,
                vector_data const& t_uv, vector_data const& br_uv,
                vector_data const& tl_uv, scalar_data const& m_temp,
                scalar_data const& r_temp, scalar_data const& t_temp)
            -> vector_partition
            {
                uint size_x = m_fg.size_x();
                uint size_y = m_fg.size_y();

                //left
                for (auto& idx_pair : boundary[0])
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    m_fg(i, j).first = m_uv(i, j).first;
                }

                //bottom
                for (auto& idx_pair : boundary[2])
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    m_fg(i, j).second = m_uv(i, j).second;
                }

                //obstacle
                for (auto& idx_pair : obstacle)
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    auto& type = flag_data[j * size_x + i];

                    if (type.test(has_fluid_east))
                        m_fg(i, j).first = m_uv(i, j).first;

                    if (type.test(has_fluid_north))
                        m_fg(i, j).second = m_uv(i, j).second;
                }

                for (auto& idx_pair : fluid)
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    compute_fg_for_cell(
                        m_fg(i, j),
                        m_uv(i, j),
                        get_left_neighbor(m_uv, l_uv, i, j),
                        get_right_neighbor(m_uv, r_uv, i, j),
                        get_bottom_neighbor(m_uv, b_uv, i, j),
                        get_top_neighbor(m_uv, t_uv, i, j),
                        get_bottomright_neighbor(m_uv, b_uv, r_uv, br_uv, i, j),
                        get_topleft_neighbor(m_uv, t_uv, l_uv, tl_uv, i, j),
                        m_temp(i, j),
                        get_right_neighbor(m_temp, r_temp, i, j),
                        get_top_neighbor(m_temp, t_temp, i, j),
                        flag_data[j * size_x + i],
                        re, gx, gy, beta, dx, dy, dt, alpha
                        );
                }

                return middle_fg;
            }
        ),
        middle_fg.get_data(CENTER),
        middle_uv.get_data(CENTER),
        left_uv.get_data(LEFT),
        right_uv.get_data(RIGHT),
        bottom_uv.get_data(BOTTOM),
        top_uv.get_data(TOP),
        bottomright_uv.get_data(BOTTOM_RIGHT),
        topleft_uv.get_data(TOP_LEFT),
        middle_temperature.get_data(CENTER),
        right_temperature.get_data(RIGHT),
        top_temperature.get_data(TOP)
    );
}

scalar_partition custom_grain_size::compute_temperature_on_fluid_cells(
    scalar_partition const& middle_temperature,
    scalar_partition const& left_temperature,
    scalar_partition const& right_temperature,
    scalar_partition const& bottom_temperature,
    scalar_partition const& top_temperature,
    vector_partition const& middle_uv, vector_partition const& left_uv,
    vector_partition const& bottom_uv,
    std::vector<std::vector<std::pair<uint, uint> > > const& boundary,
    std::vector<std::pair<uint, uint> > const& obstacle,
    std::vector<std::pair<uint, uint> > const& fluid, RealType re, RealType pr,
    RealType dx, RealType dy, RealType dt, RealType alpha)
{
    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [middle_temperature, boundary, obstacle, fluid, re, pr, dx, dy, dt,
                alpha]
            (scalar_data const& m_temp, scalar_data const& l_temp,
                scalar_data const& r_temp, scalar_data const& b_temp,
                scalar_data const& t_temp, vector_data const& m_uv,
                vector_data const& l_uv, vector_data const& b_uv)
            -> scalar_partition
            {
                uint size_x = m_temp.size_x();
                uint size_y = m_temp.size_y();

                scalar_data next(size_x, size_y);

                for (auto& idx_pair : fluid)
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    next(i, j) = compute_temperature_for_cell(
                        m_temp(i, j),
                        get_left_neighbor(m_temp, l_temp, i, j),
                        get_right_neighbor(m_temp, r_temp, i, j),
                        get_bottom_neighbor(m_temp, b_temp, i, j),
                        get_top_neighbor(m_temp, t_temp, i, j),
                        m_uv(i, j),
                        get_left_neighbor(m_uv, l_uv, i, j),
                        get_bottom_neighbor(m_uv, b_uv, i, j),
                        re, pr, dx, dy, dt, alpha);
                }

                for (auto& idx_pair : boundary[0])
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    next(i, j) = m_temp(i, j);
                }

                for (auto& idx_pair : boundary[1])
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    next(i, j) = m_temp(i, j);
                }

                for (auto& idx_pair : boundary[2])
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    next(i, j) = m_temp(i, j);
                }

                for (auto& idx_pair : boundary[3])
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    next(i, j) = m_temp(i, j);
                }

                for (auto& idx_pair : obstacle)
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    next(i, j) = m_temp(i, j);
                }

                return scalar_partition(middle_temperature.get_id(), next);
            }
        ),
        middle_temperature.get_data(CENTER),
        left_temperature.get_data(LEFT),
        right_temperature.get_data(RIGHT),
        bottom_temperature.get_data(BOTTOM),
        top_temperature.get_data(TOP),
        middle_uv.get_data(CENTER),
        left_uv.get_data(LEFT),
        bottom_uv.get_data(BOTTOM)
    );
}

scalar_partition custom_grain_size::compute_right_hand_side_on_fluid_cells(
    scalar_partition const& middle_rhs, vector_partition const& middle_fg,
    vector_partition const& left_fg, vector_partition const& bottom_fg,
    std::vector<std::pair<uint, uint> > const& fluid,
    RealType dx, RealType dy, RealType dt)
{
    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [middle_rhs, fluid, dx, dy, dt]
            (scalar_data m_rhs, vector_data const& m_fg, vector_data const& l_fg,
                vector_data const& b_fg)
            -> scalar_partition
            {
                for (auto& idx_pair : fluid)
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    m_rhs(i, j) = compute_rhs_for_cell(
                        m_fg(i, j),
                        get_left_neighbor(m_fg, l_fg, i, j),
                        get_bottom_neighbor(m_fg, b_fg, i, j),
                        dx, dy, dt);
                }

                return middle_rhs;
            }
        ),
        middle_rhs.get_data(CENTER),
        middle_fg.get_data(CENTER),
        left_fg.get_data(LEFT),
        bottom_fg.get_data(BOTTOM)
    );
}

scalar_data custom_grain_size::set_pressure_for_boundary_and_obstacles(
    hpx::shared_future<scalar_data> middle_p, hpx::shared_future<scalar_data>  left_p,
    hpx::shared_future<scalar_data> right_p, hpx::shared_future<scalar_data>  bottom_p,
    hpx::shared_future<scalar_data> top_p,
    std::vector<std::vector<std::pair<uint, uint> > > const& boundary,
    std::vector<std::pair<uint, uint> > const& obstacle,
    std::vector<std::bitset<5> > const& flag_data)
{
    auto m_p = middle_p.get();
    auto l_p = left_p.get();
    auto r_p = right_p.get();
    auto b_p = bottom_p.get();
    auto t_p = top_p.get();
    uint size_x = m_p.size_x();
    uint size_y = m_p.size_y();

    for (auto& idx_pair : boundary[0])
    {
        uint i = idx_pair.first;
        uint j = idx_pair.second;

        m_p(i, j) = m_p(i + 1, j);
    }

    for (auto& idx_pair : boundary[1])
    {
        uint i = idx_pair.first;
        uint j = idx_pair.second;

        m_p(i, j) = m_p(i - 1, j);
    }

    for (auto& idx_pair : boundary[2])
    {
        uint i = idx_pair.first;
        uint j = idx_pair.second;

        m_p(i, j) = m_p(i, j + 1);
    }

    for (auto& idx_pair : boundary[3])
    {
        uint i = idx_pair.first;
        uint j = idx_pair.second;

        m_p(i, j) = m_p(i, j - 1);
    }

    for (auto& idx_pair : obstacle)
    {
        uint i = idx_pair.first;
        uint j = idx_pair.second;

        computation::set_pressure_for_cell(
                    m_p(i, j),
                    get_left_neighbor(m_p, l_p, i, j),
                    get_right_neighbor(m_p, r_p, i, j),
                    get_bottom_neighbor(m_p, b_p, i, j),
                    get_top_neighbor(m_p, t_p, i, j),
                    flag_data[j * size_x + i]);
    }

    return m_p;
}

scalar_data custom_grain_size::sor_cycle(
    hpx::shared_future<scalar_data> middle_p, hpx::shared_future<scalar_data>  left_p,
    hpx::shared_future<scalar_data> right_p, hpx::shared_future<scalar_data>  bottom_p,
    hpx::shared_future<scalar_data> top_p, hpx::shared_future<scalar_data>  middle_rhs,
    std::vector<std::pair<uint, uint> > const& fluid,
    RealType dx_sq, RealType dy_sq, RealType part1, RealType part2)
{
    auto m_p = middle_p.get();
    auto l_p = left_p.get();
    auto r_p = right_p.get();
    auto b_p = bottom_p.get();
    auto t_p = top_p.get();
    auto m_rhs = middle_rhs.get();

    uint size_x = m_p.size_x();
    uint size_y = m_p.size_y();

    for (auto& idx_pair : fluid)
    {
        uint i = idx_pair.first;
        uint j = idx_pair.second;

        m_p(i, j)  = computation::do_sor_cycle_for_cell(
                    m_p(i, j),
                    get_left_neighbor(m_p, l_p, i, j),
                    get_right_neighbor(m_p, r_p, i, j),
                    get_bottom_neighbor(m_p, b_p, i, j),
                    get_top_neighbor(m_p, t_p, i, j),
                    m_rhs[j * size_x + i ],
                    dx_sq, dy_sq, part1, part2);
    }

    return m_p;
}

RealType custom_grain_size::compute_residual(
    hpx::shared_future<scalar_data> middle_p, hpx::shared_future<scalar_data>  left_p,
    hpx::shared_future<scalar_data> right_p, hpx::shared_future<scalar_data>  bottom_p,
    hpx::shared_future<scalar_data> top_p, hpx::shared_future<scalar_data>  middle_rhs,
    std::vector<std::pair<uint, uint> > const& fluid,
    RealType dx, RealType dy)
{
    RealType const over_dx_sq = 1./std::pow(dx, 2);
    RealType const over_dy_sq = 1./std::pow(dy, 2);
    auto m_p = middle_p.get();
    auto l_p = left_p.get();
    auto r_p = right_p.get();
    auto b_p = bottom_p.get();
    auto t_p = top_p.get();
    auto m_rhs = middle_rhs.get();
    uint size_x = m_p.size_x();
    uint size_y = m_p.size_y();

    RealType local_residual = 0;

    for (auto& idx_pair : fluid)
    {
        uint i = idx_pair.first;
        uint j = idx_pair.second;

        local_residual += compute_residual_for_cell(
            m_p(i, j),
            get_left_neighbor(m_p, l_p, i, j),
            get_right_neighbor(m_p, r_p, i, j),
            get_bottom_neighbor(m_p, b_p, i, j),
            get_top_neighbor(m_p, t_p, i, j),
            m_rhs(i, j), over_dx_sq, over_dy_sq);
    }

    return local_residual;
}

hpx::future<std::pair<vector_partition, std::pair<RealType, RealType> > >
custom_grain_size::update_velocities(
    vector_partition const& middle_uv, scalar_partition const& middle_p,
    scalar_partition const& right_p, scalar_partition const& top_p,
    vector_partition const& middle_fg,
    std::vector<std::bitset<5> > const& flag_data,
    std::vector<std::pair<uint, uint> > const& fluid,
    RealType dx, RealType dy, RealType dt)
{
    RealType const over_dx = 1./dx;
    RealType const over_dy = 1./dy;

    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [middle_uv, flag_data, fluid, over_dx, over_dy, dt]
            (vector_data m_uv, scalar_data const& m_p, scalar_data const& r_p,
                scalar_data const& t_p, vector_data const& m_fg)
            -> std::pair<vector_partition, std::pair<RealType, RealType> >
            {
                uint size_x = m_uv.size_x();
                uint size_y = m_uv.size_y();

                auto max_uv = std::make_pair(0., 0.);

                for (auto& idx_pair : fluid)
                {
                    uint i = idx_pair.first;
                    uint j = idx_pair.second;

                    auto& middle_cell = m_uv(i, j);

                    update_velocity_for_cell(
                        middle_cell,
                        m_p(i, j),
                        get_right_neighbor(m_p, r_p, i, j),
                        get_top_neighbor(m_p, t_p, i, j),
                        m_fg(i, j),
                        flag_data[j * size_x + i],
                        over_dx, over_dy, dt);

                    max_uv.first =
                       (std::abs(middle_cell.first)
                            > max_uv.first)
                        ? std::abs(middle_cell.first)
                        : max_uv.first;

                    max_uv.second =
                       (std::abs(middle_cell.second)
                            > max_uv.second)
                        ? std::abs(middle_cell.second)
                        : max_uv.second;
                }

                return std::make_pair(middle_uv, max_uv);
            }
        ),
        middle_uv.get_data(CENTER),
        middle_p.get_data(CENTER),
        right_p.get_data(RIGHT),
        top_p.get_data(TOP),
        middle_fg.get_data(CENTER)
    );
}

hpx::future<std::tuple<scalar_partition, scalar_partition, scalar_partition> >
custom_grain_size::compute_stream_vorticity_heat(
    scalar_partition const& middle_stream,
    scalar_partition const& bottom_stream,
    scalar_partition const& middle_vorticity,
    scalar_partition const& middle_heat,
    scalar_partition const& bottom_heat, vector_partition const& middle_uv,
    vector_partition const& right_uv, vector_partition const& top_uv,
    scalar_partition const& middle_temperature,
    scalar_partition const& right_temperature,
    std::vector<std::bitset<5> > const& flag_data, uint global_i, uint global_j,
    uint i_max, uint j_max, RealType re, RealType pr, RealType dx, RealType dy)
{
    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
        [middle_stream, flag_data, global_i, global_j, i_max, j_max, re, pr,
            dx, dy]
        (scalar_data stream_center, scalar_data const& stream_bottom,
            scalar_data vorticity_center, scalar_data heat_center,
            scalar_data const& heat_bottom, vector_data const& uv_center,
            vector_data const& uv_right, vector_data const& uv_top,
            scalar_data const& temp_center, scalar_data const& temp_right)
        -> std::tuple<scalar_partition, scalar_partition, scalar_partition>
        {
            uint size_x = stream_center.size_x();
            uint size_y = stream_center.size_y();

            for (uint j = 0; j < size_y; j++)
                for (uint i = 0; i < size_x; i++)
                {
                    std::bitset<5> const& cell_type = flag_data[j*size_x + i];

                    if (in_range(0, i_max, 1, j_max,
                                    global_i + i, global_j + j))
                    {
                        if (cell_type.test(4) || global_i + i == 0)
                        {
                            stream_center(i, j)  =
                                get_bottom_neighbor(
                                    stream_center, stream_bottom, i, j)
                                + uv_center(i, j).first*dy;

                            heat_center(i, j)  =
                                get_bottom_neighbor(
                                    heat_center, heat_bottom, i, j)
                                + dy * (
                                    re * pr * uv_center(i, j).first
                                    *   (get_right_neighbor(
                                            temp_center, temp_right, i, j)
                                        + temp_center(i, j)
                                        ) / 2.
                                    -   (get_right_neighbor(
                                            temp_center, temp_right, i, j)
                                        - temp_center(i, j)
                                        ) / dx
                                );
                        }
                        else
                        {
                             stream_center(i, j)  =
                                 get_bottom_neighbor(
                                    stream_center, stream_bottom, i, j) ;

                             heat_center(i, j)  =
                                 get_bottom_neighbor(
                                    heat_center, heat_bottom, i, j) ;
                        }
                    }

                    if (in_range(1, i_max - 1, 1, j_max - 1,
                                    global_i + i, global_j + j))
                    {
                        if (cell_type.test(4))
                        {
                            vector_cell const curr_cell =
                                uv_center(i, j);
                            vorticity_center(i, j)  =
                                (get_top_neighbor(uv_center, uv_top, i, j).first
                                    - curr_cell.first
                                )/dy
                                -
                                (get_right_neighbor(uv_center,
                                        uv_right, i, j).second
                                    - curr_cell.second
                                )/dx;
                        }
                        else
                            vorticity_center(i, j)  = 0;
                    }
                }

            return std::make_tuple(
                    scalar_partition(middle_stream.get_id(),stream_center),
                    scalar_partition(middle_stream.get_id(), vorticity_center),
                    scalar_partition(middle_stream.get_id(), heat_center));
        }),
        middle_stream.get_data(CENTER),
        bottom_stream.get_data(BOTTOM),
        middle_vorticity.get_data(CENTER),
        middle_heat.get_data(CENTER),
        bottom_heat.get_data(BOTTOM),
        middle_uv.get_data(CENTER),
        right_uv.get_data(RIGHT),
        top_uv.get_data(TOP),
        middle_temperature.get_data(CENTER),
        right_temperature.get_data(RIGHT));
}
}//computation

