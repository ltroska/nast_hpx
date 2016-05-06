#include "with_for_each.hpp"

#include <hpx/parallel/algorithm.hpp>
#include <hpx/parallel/algorithms/transform_reduce.hpp>

#include "util/helpers.hpp"
#include "cell_operations.hpp"

//TODO: only use get_*_neighbor where necessary

namespace computation {

vector_partition with_for_each::set_velocity_for_boundary_and_obstacles(
    vector_partition const& middle,
    vector_partition const& left,
    vector_partition const& right,
    vector_partition const& bottom,
    vector_partition const& top,
    std::vector<std::bitset<5> > const& flag_data,
    boundary_data const& type,
    boundary_data const& u,
    boundary_data const& v)
{   
    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [middle, flag_data, type, u, v]
            (vector_data next, vector_data const& l,
            vector_data const& r, vector_data const& b, vector_data const& t)
            -> vector_partition
            {
                uint size_x = next.size_x();
                uint size_y = next.size_y();

                auto range =
                    boost::irange(static_cast<uint>(0), size_x * size_y);

                // special case for the parallel algorithms: the top right most
                // cell is set wrongly, so we store the cell and set it manually
                // at the end
                vector_cell const topright_cell =
                    next(size_x - 2, size_y - 2);

                hpx::parallel::for_each(
                    hpx::parallel::par, boost::begin(range), boost::end(range),
                    [&](uint cnt)
                    {
                        uint const i = cnt%size_x;
                        uint const j = cnt/size_x;
                        
                        auto& cell_type = flag_data[j*size_x + i];
                        
                        set_velocity_for_cell(
                            next(i, j),
                            get_left_neighbor(next, l, i, j),
                            get_right_neighbor(next, r, i, j),
                            get_bottom_neighbor(next, b, i, j),
                            get_top_neighbor(next, t, i, j),
                            cell_type,
                            type, u, v);
                        
                        // special case for cells adjacent to right boundary
                        // since u is not set in the boundary cell
                        // but one cell to the left
                        if (cell_type == std::bitset<5>("01011"))
                        {
                            auto& left_cell = next(i - 1, j);

                            switch(static_cast<int>(type.right))
                            {
                                case 1 : left_cell.first = 0; break;
                                case 2 : left_cell.first = 0; break;
                                case 3 : left_cell.first =
                                    get_left_neighbor(next, l, i - 1, j).first;
                                    break;
                                case 4 : left_cell.first = u.right; break;
                            }
                        }
                        
                        //special case for cells adjacent to top boundary
                        //since v is set in bottom cell
                        if (cell_type == std::bitset<5>("01101"))
                        {
                            auto& bottom_cell = next(i, j - 1);
                            switch(static_cast<int>(type.top))
                            {
                                case 1 : bottom_cell.second = 0; break;
                                case 2 : bottom_cell.second = 0; break;
                                case 3 : bottom_cell.second =
                                    get_bottom_neighbor(next, b, i, j - 1).second;
                                    break;
                                case 4 : bottom_cell.second = v.top; break;
                            }
                        }
                    }
                );

                // restore correct value for top right cell if this partition 
                // is the top right most partition
                if (flag_data[(size_y - 1) * size_x + size_x - 2] ==
                        std::bitset<5>("01110")
                    && flag_data[(size_y - 2) * size_x + size_x - 1] ==
                            std::bitset<5>("00111"))
                {
                    if (type.top == 1)
                        next(size_x - 2, size_y - 1).first =
                            2*u.top - topright_cell.first;

                    if (type.top == 2 || type.top == 3)
                        next(size_x - 2, size_y - 1).first =
                            topright_cell.first;

                    if (type.right == 1)
                        next(size_x - 1, size_y - 2).second =
                            -topright_cell.second;

                    if (type.right == 2 || type.right == 3)
                        next(size_x - 1, size_y - 2).second =
                            topright_cell.second;
                }
                
                return vector_partition(middle.get_id(), next);
            }
        ),
        middle.get_data(CENTER),
        left.get_data(LEFT),
        right.get_data(RIGHT),
        bottom.get_data(BOTTOM),
        top.get_data(TOP)
    );        
}

scalar_partition with_for_each::set_temperature_for_boundary_and_obstacles(
    scalar_partition const& middle,
    scalar_partition const& left,
    scalar_partition const& right,
    scalar_partition const& bottom,
    scalar_partition const& top,
    std::vector<std::bitset<5> > const& flag_data,
    boundary_data const& boundary_data_type,
    boundary_data const& temperature_boundary_data,
    uint global_i, uint global_j, RealType dx, RealType dy)
{
    //TODO: add temperature for obstacle cells
    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [middle, flag_data, boundary_data_type, temperature_boundary_data,
             global_i, global_j, dx, dy]
            (scalar_data next, scalar_data const& l,
            scalar_data const& r, scalar_data const& b, scalar_data const& t)
            -> scalar_partition
            {
                uint size_x = next.size_x();
                uint size_y = next.size_y();

                auto range =
                    boost::irange(static_cast<uint>(0), size_x * size_y);

                hpx::parallel::for_each(
                    hpx::parallel::par, boost::begin(range), boost::end(range),
                    [&](uint cnt)
                    {
                        uint const i = cnt%size_x;
                        uint const j = cnt/size_x;

                        scalar_cell& middle = next(i, j);
                        scalar_cell const left =
                            get_left_neighbor(next, l, i, j);
                        scalar_cell const right =
                            get_right_neighbor(next, r, i, j);
                        scalar_cell const bottom =
                            get_bottom_neighbor(next, b, i, j);
                        scalar_cell const top =
                            get_top_neighbor(next, t, i, j);

                        std::bitset<5> cell_type = flag_data[j*size_x + i];
                        
                        set_temperature_for_cell(middle, left, right,
                            bottom, top, boundary_data_type,
                            temperature_boundary_data, cell_type,
                            global_i + i, global_j + j, dx, dy);             
                    }
                );
                
                return scalar_partition(middle.get_id(), next);
            }
        ),
        middle.get_data(CENTER),
        left.get_data(LEFT),
        right.get_data(RIGHT),
        bottom.get_data(BOTTOM),
        top.get_data(TOP)
    );        
}


vector_partition with_for_each::compute_fg_on_fluid_cells(
    vector_partition const& middle_uv,
    vector_partition const& left_uv, vector_partition const& right_uv,
    vector_partition const& bottom_uv, vector_partition const& top_uv,
    vector_partition const& bottomright_uv,
    vector_partition const& topleft_uv,
    scalar_partition const& middle_temperature,
    scalar_partition const& right_temperature,
    scalar_partition const& top_temperature,
    std::vector<std::bitset<5> > const& flag_data,
    RealType re,
    RealType gx, RealType gy, RealType beta,
    RealType dx, RealType dy, RealType dt, RealType alpha)
{           
    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [middle_uv, flag_data, re, gx, gy, beta, dx, dy, dt, alpha]
            (vector_data const& m_uv, vector_data const& l_uv,
                vector_data const& r_uv, vector_data const& b_uv,
                vector_data const& t_uv, vector_data const& br_uv,
                vector_data const& tl_uv, scalar_data const& m_temp,
                scalar_data const& r_temp, scalar_data const& t_temp)
            -> vector_partition
            {
                uint size_x = m_uv.size_x();
                uint size_y = m_uv.size_y();
           
                vector_data next(size_x, size_y);
                
                auto range =
                    boost::irange(static_cast<uint>(0), size_x * size_y);
                
                hpx::parallel::for_each(
                    hpx::parallel::par, boost::begin(range), boost::end(range),
                [&](uint cnt)
                {
                    uint const i = cnt%size_x;
                    uint const j = cnt/size_x;

                    compute_fg_for_cell(
                        next(i, j), m_uv(i, j),
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
                        re, gx, gy, beta, dx, dy, dt, alpha);
                   }
                );                  
                
                return vector_partition(middle_uv.get_id(), next);
            }
        ),
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

scalar_partition with_for_each::compute_temperature_on_fluid_cells(
    scalar_partition const& middle_temperature,
    scalar_partition const& left_temperature,
    scalar_partition const& right_temperature,
    scalar_partition const& bottom_temperature,
    scalar_partition const& top_temperature,
    vector_partition const& middle_uv,
    vector_partition const& left_uv,
    vector_partition const& bottom_uv,
    std::vector<std::bitset<5> > const& flag_data,
    RealType re, RealType pr, RealType dx, RealType dy, RealType dt,
    RealType alpha        
)
{                      
    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [middle_temperature, flag_data, re, pr, dx, dy, dt, alpha]
            ( scalar_data const& m_temp, scalar_data const& l_temp,
                scalar_data const& r_temp, scalar_data const& b_temp,
                scalar_data const& t_temp, vector_data const& m_uv,
                vector_data const& l_uv, vector_data const& b_uv)
            -> scalar_partition
            {
                uint size_x = m_temp.size_x();
                uint size_y = m_temp.size_y();
                   
                scalar_data next(m_temp);

                auto range =
                    boost::irange(static_cast<uint>(0), size_x * size_y);
                
                hpx::parallel::for_each(
                    hpx::parallel::par, boost::begin(range), boost::end(range),
                    [&](uint cnt)
                    {
                        uint const i = cnt%size_x;
                        uint const j = cnt/size_x;
                        
                        compute_temperature_for_cell(
                            next(i, j),
                            m_temp(i, j),
                            get_left_neighbor(m_temp, l_temp, i, i),
                            get_right_neighbor(m_temp, r_temp, i, j),
                            get_bottom_neighbor(m_temp, b_temp, i, j),
                            get_top_neighbor(m_temp, t_temp, i, j),
                            m_uv(i, j),
                            get_left_neighbor(m_uv, l_uv, i, j),
                            get_bottom_neighbor(m_uv, b_uv, i, j),
                            flag_data[j * size_x + i],
                            re, pr, dx, dy, dt, alpha);
                    }
                );

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



scalar_partition with_for_each::compute_right_hand_side_on_fluid_cells(
    vector_partition const& middle_fg, vector_partition const& left_fg,
    vector_partition const& bottom_fg,
    std::vector<std::bitset<5> > const& flag_data,
    RealType dx, RealType dy, RealType dt)
{  
    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [middle_fg, flag_data, dx, dy, dt]
            (vector_data const& m_fg, vector_data const& l_fg,
                vector_data const& b_fg)
            -> scalar_partition
            {
                uint size_x = m_fg.size_x();
                uint size_y = m_fg.size_y();
                
                scalar_data next(size_x, size_y);
                      
                auto range =
                    boost::irange(static_cast<uint>(0), size_x * size_y);
                
                hpx::parallel::for_each(
                    hpx::parallel::par, boost::begin(range), boost::end(range),
                    [&](uint cnt)
                    {
                        uint const i = cnt%size_x;
                        uint const j = cnt/size_x;
                        
                        compute_rhs_for_cell(
                            next(i, j),
                            m_fg(i, j),
                            get_left_neighbor(m_fg, l_fg, i, j),
                            get_bottom_neighbor(m_fg, b_fg, i, j),
                            flag_data[j * size_x + i],
                            dx, dy, dt);
                    }
                );                
                
                return scalar_partition(middle_fg.get_id(), next);
            }
        ),
        middle_fg.get_data(CENTER),
        left_fg.get_data(LEFT),
        bottom_fg.get_data(BOTTOM)
    );    
}

scalar_partition with_for_each::set_pressure_for_boundary_and_obstacles(
        scalar_partition const& middle_p, scalar_partition const& left_p,
        scalar_partition const& right_p, scalar_partition const& bottom_p,
        scalar_partition const& top_p,
        std::vector<std::bitset<5> > const& flag_data)
{   
   return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [middle_p, flag_data]
            (scalar_data next, scalar_data const& l_p,
                scalar_data const& r_p, scalar_data const& b_p,
                scalar_data const& t_p)
            -> scalar_partition
            {
                uint size_x = next.size_x();
                uint size_y = next.size_y();        
                        
                auto range =
                    boost::irange(static_cast<uint>(0), size_x * size_y);
                
                hpx::parallel::for_each(
                    hpx::parallel::par, boost::begin(range), boost::end(range),
                    [&](uint cnt)
                    {
                        uint const i = cnt%size_x;
                        uint const j = cnt/size_x;

                        set_pressure_for_cell(
                            next(i, j),
                            get_left_neighbor(next, l_p, i, j),
                            get_right_neighbor(next, r_p, i, j),
                            get_bottom_neighbor(next, b_p, i, j),
                            get_top_neighbor(next, t_p, i, j),
                            flag_data[j * size_x + i]);

                    }
                );                                            
                
                return scalar_partition(middle_p.get_id(), next);
            }
        ),
        middle_p.get_data(CENTER),
        left_p.get_data(LEFT),
        right_p.get_data(RIGHT),
        bottom_p.get_data(BOTTOM),
        top_p.get_data(TOP)
    );
}

scalar_partition with_for_each::sor_cycle(scalar_partition const& middle_p,
        scalar_partition const& left_p, scalar_partition const& right_p,
        scalar_partition const& bottom_p, scalar_partition const& top_p,
        scalar_partition const& middle_rhs, 
        std::vector<std::bitset<5> > const& flag_data,
        RealType omega, RealType dx, RealType dy)
{
    RealType const dx_sq = std::pow(dx, 2);
    RealType const dy_sq = std::pow(dy, 2);
    RealType const part1 = 1. - omega;
    RealType const part2 = omega * dx_sq * dy_sq / (2. * (dx_sq + dy_sq));

    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [middle_p, flag_data, dx_sq, dy_sq, part1, part2]
            (scalar_data m_p, scalar_data const& l_p, scalar_data const& r_p,
                scalar_data const& b_p, scalar_data const& t_p,
                scalar_data const& m_rhs)
            -> scalar_partition
            {
                uint size_x = m_p.size_x();
                uint size_y = m_p.size_y();

                auto range =
                    boost::irange(static_cast<uint>(0), size_x * size_y);
                
                hpx::parallel::for_each(
                    hpx::parallel::par, boost::begin(range), boost::end(range),
                    [&](uint cnt)
                    {
                        uint const i = cnt%size_x;
                        uint const j = cnt/size_x;

                        do_sor_cycle_for_cell(
                            m_p(i, j),
                            get_left_neighbor(m_p, l_p, i, j),
                            get_right_neighbor(m_p, r_p, i, j),
                            get_bottom_neighbor(m_p, b_p, i, j),
                            get_top_neighbor(m_p, t_p, i, j),
                            m_rhs(i, j),
                            flag_data[j * size_x + i],
                            dx_sq, dy_sq, part1, part2);
                    }
                );    
                
                return scalar_partition(middle_p.get_id(), m_p);
            }
        ),
        middle_p.get_data(CENTER),
        left_p.get_data(LEFT),
        right_p.get_data(RIGHT),
        bottom_p.get_data(BOTTOM),
        top_p.get_data(TOP),
        middle_rhs.get_data(CENTER)
    );    
}

hpx::future<RealType> with_for_each::compute_residual(
    scalar_partition const& middle_p, scalar_partition const& left_p,
    scalar_partition const& right_p, scalar_partition const& bottom_p,
    scalar_partition const& top_p, scalar_partition const& middle_rhs,
    std::vector<std::bitset<5> > const& flag_data, RealType dx, RealType dy)
{
    RealType const over_dx_sq = 1./std::pow(dx, 2);
    RealType const over_dy_sq = 1./std::pow(dy, 2);
   
    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [flag_data, over_dx_sq, over_dy_sq]
            (scalar_data const& m_p, scalar_data const& l_p,
                scalar_data const& r_p, scalar_data const& b_p,
                scalar_data const& t_p, scalar_data const& m_rhs)
            -> RealType
            {
                uint size_x = m_p.size_x();
                uint size_y = m_p.size_y();
                                
                auto range =
                    boost::irange(static_cast<uint>(0), size_x * size_y);
                                
                return hpx::parallel::transform_reduce(
                    hpx::parallel::par, boost::begin(range), boost::end(range),
                    [&](uint cnt)
                        -> RealType
                    {
                        uint const i = cnt%size_x;
                        uint const j = cnt/size_x;

                        return compute_residual_for_cell(
                            m_p(i, j),
                            get_left_neighbor(m_p, l_p, i, j),
                            get_right_neighbor(m_p, r_p, i, j),
                            get_bottom_neighbor(m_p, b_p, i, j),
                            get_top_neighbor(m_p, t_p, i, j),
                            m_rhs(i, j),
                            flag_data[j * size_x + i],
                            over_dx_sq, over_dy_sq);
                    },
                    0.,
                    [](RealType a, RealType b) -> RealType {return a + b;}
                    );                
            }
        ),
        middle_p.get_data(CENTER),
        left_p.get_data(LEFT),
        right_p.get_data(RIGHT),
        bottom_p.get_data(BOTTOM),
        top_p.get_data(TOP),
        middle_rhs.get_data(CENTER)
    );   
}

hpx::future<std::pair<vector_partition, std::pair<RealType, RealType> > > 
with_for_each::update_velocities(
    vector_partition const& middle_uv, scalar_partition const& middle_p,
    scalar_partition const& right_p, scalar_partition const& top_p, 
    vector_partition const& middle_fg,
    std::vector<std::bitset<5> > const& flag_data, RealType dx, RealType dy,
    RealType dt)
{
    RealType const over_dx = 1./dx;
    RealType const over_dy = 1./dy;
       
    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [middle_uv, flag_data, over_dx, over_dy, dt]
            (vector_data next, scalar_data const& m_p, scalar_data const& r_p,
                scalar_data const& t_p, vector_data const& m_fg)
            -> std::pair<vector_partition, std::pair<RealType, RealType> >
            {
                uint size_x = next.size_x();
                uint size_y = next.size_y();
                                
                auto range =
                    boost::irange(static_cast<uint>(0), size_x * size_y);   
                
                // we transform the velocity cells first ( = update them to the
                // new values) and then reduce to retrieve the maximal
                // velocities
                vector_cell max_uv = hpx::parallel::transform_reduce(
                    hpx::parallel::par, boost::begin(range), boost::end(range),
                    [&](uint cnt)
                        -> vector_cell
                    {
                        uint const i = cnt%size_x;
                        uint const j = cnt/size_x;

                        vector_cell& middle_uv = next(i, j);
                        
                        update_velocity_for_cell(
                            middle_uv,
                            m_p(i, j),
                            get_right_neighbor(m_p, r_p, i, j),
                            get_top_neighbor(m_p, t_p, i, j),
                            m_fg(i, j),
                            flag_data[j * size_x + i],
                            over_dx, over_dy, dt);
                        
                        return middle_uv;
                    },
                    vector_cell(0, 0),
                    [](vector_cell a, vector_cell b) -> vector_cell
                    {
                        return vector_cell(
                            std::abs(a.first) > std::abs(b.first)
                                ? std::abs(a.first) : std::abs(b.first),
                            std::abs(a.second) > std::abs(b.second)
                                ? std::abs(a.second) : std::abs(b.second));
                    }
                    );
                    
                    return std::make_pair(
                        vector_partition(middle_uv.get_id(), next),
                        std::make_pair(max_uv.first, max_uv.second));
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
with_for_each::compute_stream_vorticity_heat(
    scalar_partition const& middle_stream,
    scalar_partition const& bottom_stream,
    scalar_partition const& middle_vorticity,
    scalar_partition const& middle_heat,
    scalar_partition const& bottom_heat, vector_partition const& middle_uv,
    vector_partition const& right_uv, vector_partition const& top_uv,
    scalar_partition const& middle_temperature,
    scalar_partition const& right_temperature,
    std::vector<std::bitset<5> > const& flag_data, 
    uint global_i, uint global_j, uint i_max, uint j_max, RealType re,
    RealType pr, RealType dx, RealType dy)
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

            //TODO: use parallel algorithm here            
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
                                get_bottom_neighbor(stream_center,
                                    stream_bottom, i, j) 
                                + uv_center(i, j).first*dy;
                            
                            heat_center(i, j)  =
                                get_bottom_neighbor(heat_center,
                                    heat_bottom, i, j) 
                                + dy * (re * pr * uv_center(i, j).first
                                * (get_right_neighbor(temp_center,
                                                        temp_right, i, j) 
                                    + temp_center(i, j) ) / 2.
                                - (get_right_neighbor(temp_center,
                                        temp_right, i, j) 
                                    - temp_center(i, j) ) / dx);
                        }
                        else
                        {
                             stream_center(i, j)  =
                                get_bottom_neighbor(stream_center,
                                    stream_bottom, i, j) ;
                             
                             heat_center(i, j)  =
                                get_bottom_neighbor(heat_center,
                                    heat_bottom, i, j) ;
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
                                    - curr_cell.first)/dy
                                - (get_right_neighbor(uv_center,
                                            uv_right, i, j).second 
                                    - curr_cell.second)/dx;
                        }
                        else
                            vorticity_center(i, j)  = 0;
                    }
                }
            
            return std::make_tuple(
                    scalar_partition(middle_stream.get_id(), stream_center),
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

