#include "custom_grain_size.hpp"

#include "util/helpers.hpp"
#include "stencils.hpp"

//TODO: only use get_*_neighbor where necessary

namespace computation {

vector_partition custom_grain_size::set_velocity_for_boundary_and_obstacles(
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
    //do local computation first
    hpx::future<vector_data> next_middle = 
        hpx::dataflow(
            hpx::launch::async,
            hpx::util::unwrapped(
                [middle, flag_data, type, u, v]
                (vector_data m) -> vector_data
                {
                    uint size_x = m.size_x();
                    uint size_y = m.size_y();

                    for (uint j = 1; j < size_y - 1; j++)
                        for (uint i = 1; i < size_x - 1; i++)
                            set_velocity_for_cell(
                                m.get_cell_ref(i, j),
                                m.get_cell(i - 1, j),
                                m.get_cell(i + 1, j),
                                m.get_cell(i, j - 1),
                                m.get_cell(i, j + 1),
                                flag_data[j * size_x + i],
                                type, u, v);

                    return m;
                }
            ),
            middle.get_data(CENTER)
        );

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
                                
                //left and right
                for (uint j = 0; j < size_y; j++)
                {
                    uint i = 0;
                    set_velocity_for_cell(
                        next.get_cell_ref(i, j),
                        get_left_neighbor(next, l, i, j),
                        next.get_cell(i + 1, j),
                        get_bottom_neighbor(next, b, i, j),
                        get_top_neighbor(next, t, i, j),
                        flag_data[j * size_x + i],
                        type, u, v);
                              
                    i = size_x - 1;
                    auto& cell_type = flag_data[j * size_x + i];
                    
                    set_velocity_for_cell(
                        next.get_cell_ref(i, j),
                        next.get_cell(i - 1, j),
                        get_right_neighbor(next, r, i, j),
                        get_bottom_neighbor(next, b, i, j),
                        get_top_neighbor(next, t, i, j),
                        cell_type,
                        type, u, v);
                    
                    //special case for cells adjacent to right boundary
                    ///since u is set in left cell
                    if (cell_type == std::bitset<5>("01011"))
                    {
                        auto& left_cell = next.get_cell_ref(i - 1, j);
                        
                        switch((int)type.right)
                        {
                            case 1 : left_cell.first = 0; break;
                            case 2 : left_cell.first = 0; break;
                            case 3 : left_cell.first =
                                get_left_neighbor(next, l, i - 1, j).first;
                                break;
                            case 4 : left_cell.first = u.right; break;
                        }
                    }
                }    
                
                //bottom and top
                for (uint i = 0; i < size_x; i++)
                {
                    uint j = 0;
                    set_velocity_for_cell(
                        next.get_cell_ref(i, j),
                        get_left_neighbor(next, l, i, j),
                        get_right_neighbor(next, r, i, j),
                        get_bottom_neighbor(next, b, i, j),
                        next.get_cell(i, j + 1),
                        flag_data[j * size_x + i],
                        type, u, v);
                     
                    
                    j = size_y - 1;
                    auto& cell_type = flag_data[j * size_x + i];
                    
                    set_velocity_for_cell(
                        next.get_cell_ref(i, j),
                        get_left_neighbor(next, l, i, j),
                        get_right_neighbor(next, r, i, j),
                        next.get_cell(i, j - 1),
                        get_top_neighbor(next, t, i, j),
                        cell_type,
                        type, u, v);
                    
                    //special case for cells adjacent to top boundary
                    //since v is set in bottom cell
                    if (cell_type == std::bitset<5>("01101"))
                    {
                        auto& bottom_cell = next.get_cell_ref(i, j - 1);
                        switch((int)type.top)
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
                
                return vector_partition(middle.get_id(), next);
            }
        ),
        std::move(next_middle),
        left.get_data(LEFT),
        right.get_data(RIGHT),
        bottom.get_data(BOTTOM),
        top.get_data(TOP)
    );        
}

void custom_grain_size::set_velocity_for_cell(vector_cell& middle,
    vector_cell const& left, vector_cell const& right,
    vector_cell const& bottom, vector_cell const& top,
    std::bitset<5> const& cell_type,
    boundary_data const& type,
    boundary_data const& u,
    boundary_data const& v)
{
    //obstacles
    //north
    if (cell_type == std::bitset<5>("00001"))
    {
        middle.second = 0;
        middle.first = -top.first;
    }

    //east
    else if (cell_type == std::bitset<5>("01000"))
    {
        middle.first = 0;
        middle.second = -right.second;
    }

    //south
    else if (cell_type == std::bitset<5>("00010"))
    {
        middle.first = -bottom.first;
    }

    else if (cell_type.test(4) && !cell_type.test(0))
    {
        middle.second = 0;          
    }

    //west
    else if (cell_type == std::bitset<5>("00100"))
    {
        middle.second = -left.second;
    }

    else if (cell_type.test(4) && !cell_type.test(3))
    {
            middle.first = 0;
    }

    //north east
    else if (cell_type == std::bitset<5>("01001"))
    {
        middle.first = 0;
        middle.second = 0;
    }

    //south east
    else if (cell_type == std::bitset<5>("01010"))
    {
        middle.first = 0;
        middle.second = -right.second;
    }

    //south west
    else if (cell_type == std::bitset<5>("00110"))
    {
        middle.first = -bottom.first;
        middle.second = -left.second;
    }

    //north west
    else if (cell_type == std::bitset<5>("00101"))
    {
        middle.first = -top.first;
        middle.second = 0;
    }
    
    //boundary cells
    //left
    else if (cell_type == std::bitset<5>("00111"))
        switch((int)type.left)
        {
        case 1: middle.first = 0;
                middle.second = -right.second;
                break;
        case 2: middle.first = 0;
                middle.second = right.second;
                break;
        case 3: middle.first = right.first;
                middle.second = right.second;
                break;
        case 4: middle.first = u.left;
                middle.second = 2*v.left - right.second;
                break;  
        }
    
    //right
    else if (cell_type == std::bitset<5>("01011"))
        switch((int)type.right)
        {
        case 1: middle.second = -left.second;
                break;
        case 2: middle.second = left.second;
                break;
        case 3: middle.second = left.second;
                break;
        case 4: middle.second = 2*v.right - left.second;
                break;  
        }      

    //bottom
    else if (cell_type == std::bitset<5>("01110"))
        switch((int)type.bottom)
        {
        case 1: middle.first = -top.first;
                middle.second = 0;
                break;
        case 2: middle.first = top.first;
                middle.second = 0;
                break;
        case 3: middle.first = top.first;
                middle.second = top.second;
                break;
        case 4: middle.first = 2*u.bottom - top.first;
                middle.second = v.bottom;
                break;  
        }  

    //top
    else if (cell_type == std::bitset<5>("01101"))
        switch((int)type.top)
        {
        case 1: middle.first = 2*u.top - bottom.first;
                break;
        case 2: middle.first = bottom.first;
                break;
        case 3: middle.first = bottom.first;
                break;
        case 4: middle.first = 2*u.top - bottom.first;
                break;  
        } 
}  

scalar_partition custom_grain_size::set_temperature_for_boundary_and_obstacles(
        scalar_partition const& middle,
        scalar_partition const& left,
        scalar_partition const& right,
        scalar_partition const& bottom,
        scalar_partition const& top,
        std::vector<std::bitset<5> > const& flag_data,
        boundary_data const& boundary_data_type,
        boundary_data const& temperature_boundary_data,
        RealType dx, RealType dy)
{
    //TODO: add temperature for obstacle cells
    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [middle, flag_data, boundary_data_type, temperature_boundary_data,
             dx, dy]
            (scalar_data next, scalar_data const& l,
            scalar_data const& r, scalar_data const& b, scalar_data const& t)
            -> scalar_partition
            {
                uint size_x = next.size_x();
                uint size_y = next.size_y();
             
                //left and right
                for (uint j = 0; j < size_y; j++)
                {
                    uint i = 0;
                    set_temperature_for_cell(
                        next.get_cell_ref(i, j),
                        get_left_neighbor(next, l, i, j),
                        get_right_neighbor(next, r, i, j),
                        get_bottom_neighbor(next, b, i, j),
                        get_top_neighbor(next, t, i, j),
                        boundary_data_type, temperature_boundary_data,
                        flag_data[j * size_x + i], i, j,dx, dy);
                     
                    
                    i = size_x - 1;
                    set_temperature_for_cell(
                        next.get_cell_ref(i, j),
                        get_left_neighbor(next, l, i, j),
                        get_right_neighbor(next, r, i, j),
                        get_bottom_neighbor(next, b, i, j),
                        get_top_neighbor(next, t, i, j),
                        boundary_data_type, temperature_boundary_data,
                        flag_data[j * size_x + i],i, j, dx, dy);
                }    
                
                //bottom and top
                for (uint i = 0; i < size_x; i++)
                {
                    uint j = 0;
                    set_temperature_for_cell(
                        next.get_cell_ref(i, j),
                        get_left_neighbor(next, l, i, j),
                        get_right_neighbor(next, r, i, j),
                        get_bottom_neighbor(next, b, i, j),
                        get_top_neighbor(next, t, i, j),
                        boundary_data_type, temperature_boundary_data,
                        flag_data[j * size_x + i], i, j,dx, dy);
                     
                    
                    j = size_y - 1;
                    set_temperature_for_cell(
                        next.get_cell_ref(i, j),
                        get_left_neighbor(next, l, i, j),
                        get_right_neighbor(next, r, i, j),
                        get_bottom_neighbor(next, b, i, j),
                        get_top_neighbor(next, t, i, j),
                        boundary_data_type, temperature_boundary_data,
                        flag_data[j * size_x + i], i, j,dx, dy);
                }                     
                
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

void custom_grain_size::set_temperature_for_cell(scalar_cell& middle,
        scalar_cell const& left, scalar_cell const& right,
        scalar_cell const& bottom, scalar_cell const& top,
        boundary_data const& boundary_data_type,
        boundary_data const& temperature_boundary_data,
        std::bitset<5> const& cell_type,
        uint i, uint j,
        RealType dx, RealType dy
        )
{
    //left
    if (cell_type == std::bitset<5>("00111"))
        switch((int)boundary_data_type.left)
        {
        case 1: middle.value = 2*temperature_boundary_data.left - right.value; break;            
        case 2: middle.value = right.value + dx*temperature_boundary_data.left*((j-0.5)*dy); break;
        }
    
    //right
    if (cell_type == std::bitset<5>("01011"))
        switch((int)boundary_data_type.right)
        {
        case 1: middle.value = 2*temperature_boundary_data.right - left.value; break;            
        case 2: middle.value = left.value + dx*temperature_boundary_data.right*((j-0.5)*dy); break;
        }

    //bottom
    if (cell_type == std::bitset<5>("01110"))
        switch((int)boundary_data_type.bottom)
        {
        case 1: middle.value = 2*temperature_boundary_data.bottom - top.value; break;            
        case 2: middle.value = top.value + dy*temperature_boundary_data.bottom*((i-0.5)*dx); break;
        }

    //top
    if (cell_type == std::bitset<5>("01101"))
        switch((int)boundary_data_type.top)
        {
        case 1: middle.value = 2*temperature_boundary_data.top - bottom.value; break;            
        case 2: middle.value = bottom.value + dy*temperature_boundary_data.top*((i-0.5)*dx); break;
        }    
}

vector_partition custom_grain_size::compute_fg_on_fluid_cells(
      vector_partition const& middle_uv,
      vector_partition const& left_uv, vector_partition const& right_uv,
      vector_partition const& bottom_uv, vector_partition const& top_uv,
      vector_partition const& bottomright_uv, vector_partition const& topleft_uv,
      scalar_partition const& middle_temperature,
      scalar_partition const& right_temperature,
      scalar_partition const& top_temperature,
      std::vector<std::bitset<5> > const& flag_data,
      RealType re,
      RealType gx, RealType gy, RealType beta,
      RealType dx, RealType dy, RealType dt, RealType alpha)
{           
    
    hpx::shared_future<scalar_data> middle_temperature_data =
        middle_temperature.get_data(CENTER);
    
    hpx::shared_future<vector_data> middle_uv_data = middle_uv.get_data(CENTER);
    
    //do local computation first
     hpx::future<vector_data> next_fg_middle = 
        hpx::dataflow(
            hpx::launch::async,
            hpx::util::unwrapped(
                [middle_uv, flag_data, re, gx, gy, beta, dx, dy, dt, alpha]
                (vector_data const& m_uv,
                    scalar_data const& m_temp)
                -> vector_data
                {
                    uint size_x = m_uv.size_x();
                    uint size_y = m_uv.size_y();

                    vector_data next(size_x, size_y);

                    for (uint j = 1; j < size_y - 1; j++)
                        for (uint i = 1; i < size_x - 1; i++)
                            compute_fg_for_cell(
                                    next.get_cell_ref(i, j),
                                    m_uv.get_cell(i, j),
                                    m_uv.get_cell(i - 1, j),
                                    m_uv.get_cell(i + 1, j),
                                    m_uv.get_cell(i, j - 1),
                                    m_uv.get_cell(i, j + 1),
                                    m_uv.get_cell(i + 1, j - 1),
                                    m_uv.get_cell(i - 1, j + 1),
                                    m_temp.get_cell(i, j),
                                    m_temp.get_cell(i + 1, j),
                                    m_temp.get_cell(i, j + 1),
                                    flag_data[j * size_x + i],
                                    re, gx, gy, beta, dx, dy, dt, alpha);

                    return next;
                }
            ),
            middle_uv_data,
            middle_temperature_data
        );

    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [middle_uv, flag_data, re, gx, gy, beta, dx, dy, dt, alpha]
            (vector_data next, vector_data const& m_uv, vector_data const& l_uv,
            vector_data const& r_uv, vector_data const& b_uv, vector_data const& t_uv,
            vector_data const& br_uv, vector_data const& tl_uv,
            scalar_data const& m_temp, scalar_data const& r_temp,
            scalar_data const& t_temp)
            -> vector_partition
            {
                uint size_x = next.size_x();
                uint size_y = next.size_y();
                                
                //left and right
                for (uint j = 0; j < size_y; j++)
                {
                    uint i = 0;                     
                    compute_fg_for_cell(
                        next.get_cell_ref(i, j),
                        m_uv.get_cell(i, j),
                        get_left_neighbor(m_uv, l_uv, i, j),
                        get_right_neighbor(m_uv, r_uv, i, j),
                        get_bottom_neighbor(m_uv, b_uv, i, j),
                        get_top_neighbor(m_uv, t_uv, i, j),
                        get_bottomright_neighbor(m_uv, b_uv, r_uv, br_uv, i, j),
                        get_topleft_neighbor(m_uv, t_uv, l_uv, tl_uv, i, j),
                        m_temp.get_cell(i, j),
                        get_right_neighbor(m_temp, r_temp, i, j),
                        get_top_neighbor(m_temp, t_temp, i, j),
                        flag_data[j * size_x + i],
                        re, gx, gy, beta, dx, dy, dt, alpha);
                     
                    
                    i = size_x - 1;
                    compute_fg_for_cell(
                        next.get_cell_ref(i, j),
                        m_uv.get_cell(i, j),
                        get_left_neighbor(m_uv, l_uv, i, j),
                        get_right_neighbor(m_uv, r_uv, i, j),
                        get_bottom_neighbor(m_uv, b_uv, i, j),
                        get_top_neighbor(m_uv, t_uv, i, j),
                        get_bottomright_neighbor(m_uv, b_uv, r_uv, br_uv, i, j),
                        get_topleft_neighbor(m_uv, t_uv, l_uv, tl_uv, i, j),
                        m_temp.get_cell(i, j),
                        get_right_neighbor(m_temp, r_temp, i, j),
                        get_top_neighbor(m_temp, t_temp, i, j),
                        flag_data[j * size_x + i],
                        re, gx, gy, beta, dx, dy, dt, alpha);
                }    
                
                //bottom and top
                for (uint i = 0; i < size_x; i++)
                {
                    uint j = 0;
                    compute_fg_for_cell(
                        next.get_cell_ref(i, j),
                        m_uv.get_cell(i, j),
                        get_left_neighbor(m_uv, l_uv, i, j),
                        get_right_neighbor(m_uv, r_uv, i, j),
                        get_bottom_neighbor(m_uv, b_uv, i, j),
                        get_top_neighbor(m_uv, t_uv, i, j),
                        get_bottomright_neighbor(m_uv, b_uv, r_uv, br_uv, i, j),
                        get_topleft_neighbor(m_uv, t_uv, l_uv, tl_uv, i, j),
                        m_temp.get_cell(i, j),
                        get_right_neighbor(m_temp, r_temp, i, j),
                        get_top_neighbor(m_temp, t_temp, i, j),
                        flag_data[j * size_x + i],
                        re, gx, gy, beta, dx, dy, dt, alpha);
                     
                    
                    j = size_y - 1;
                    compute_fg_for_cell(
                        next.get_cell_ref(i, j),
                        m_uv.get_cell(i, j),
                        get_left_neighbor(m_uv, l_uv, i, j),
                        get_right_neighbor(m_uv, r_uv, i, j),
                        get_bottom_neighbor(m_uv, b_uv, i, j),
                        get_top_neighbor(m_uv, t_uv, i, j),
                        get_bottomright_neighbor(m_uv, b_uv, r_uv, br_uv, i, j),
                        get_topleft_neighbor(m_uv, t_uv, l_uv, tl_uv, i, j),
                        m_temp.get_cell(i, j),
                        get_right_neighbor(m_temp, r_temp, i, j),
                        get_top_neighbor(m_temp, t_temp, i, j),
                        flag_data[j * size_x + i],
                        re, gx, gy, beta, dx, dy, dt, alpha);
                }                     
                
                return vector_partition(middle_uv.get_id(), next);
            }
        ),
        std::move(next_fg_middle),
        middle_uv_data,
        left_uv.get_data(LEFT),
        right_uv.get_data(RIGHT),
        bottom_uv.get_data(BOTTOM),
        top_uv.get_data(TOP),
        bottomright_uv.get_data(BOTTOM_RIGHT),
        topleft_uv.get_data(TOP_LEFT),
        middle_temperature_data,
        right_temperature.get_data(RIGHT),
        top_temperature.get_data(TOP)
    );        
}

void custom_grain_size::compute_fg_for_cell(vector_cell& middle_fg,
        vector_cell const& middle_uv, vector_cell const& left_uv,
        vector_cell const& right_uv, vector_cell const& bottom_uv,
        vector_cell const& top_uv, vector_cell const& bottomright_uv,
        vector_cell const& topleft_uv, scalar_cell const& middle_temperature,
        scalar_cell const& right_temperature, scalar_cell const& top_temperature,
        std::bitset<5> const& type,
        RealType re, RealType gx, RealType gy, RealType beta,
        RealType dx, RealType dy, RealType dt, RealType alpha)
{   
    //north
    if (!type.test(4) && type.test(0))
        middle_fg.second = middle_uv.second;

    //east
    if (!type.test(4) && type.test(3))
        middle_fg.first = middle_uv.first;

    if (type.test(4))
    {
        if (!type.test(3))
            middle_fg.first = middle_uv.first;
        else
        {
            middle_fg.first = middle_uv.first + dt * (
                            1./re * (second_derivative_fwd_bkwd_x(right_uv.first, middle_uv.first, left_uv.first, dx)
                                        + second_derivative_fwd_bkwd_y(top_uv.first, middle_uv.first, bottom_uv.first, dy))

                            - first_derivative_of_square_x(right_uv.first, middle_uv.first, left_uv.first, dx, alpha)
                            - first_derivative_of_product_y(right_uv.second, middle_uv.second, bottom_uv.second, bottomright_uv.second,
                                                            bottom_uv.first, middle_uv.first, top_uv.first, dy, alpha)

                            + gx
                            )
                            - beta * dt / 2. * (middle_temperature.value + right_temperature.value) * gx;
        }

        if (!type.test(0))
            middle_fg.second = middle_uv.second;
        else
        {
            middle_fg.second = middle_uv.second + dt * (
                            1./re * (second_derivative_fwd_bkwd_x(right_uv.second, middle_uv.second, left_uv.second, dx)
                                        + second_derivative_fwd_bkwd_y(top_uv.second, middle_uv.second, bottom_uv.second, dy))

                            - first_derivative_of_product_x(left_uv.first, middle_uv.first, top_uv.first, topleft_uv.first,
                                                            left_uv.second, middle_uv.second, right_uv.second, dx, alpha)
                            - first_derivative_of_square_y(top_uv.second, middle_uv.second, bottom_uv.second, dy, alpha)

                            + gy
                            )
                            - beta * dt / 2. * (middle_temperature.value + top_temperature.value) * gy;
        }
    }
}

scalar_partition custom_grain_size::compute_temperature_on_fluid_cells(
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
    hpx::shared_future<scalar_data> middle_temperature_data =
        middle_temperature.get_data(CENTER);
    
    hpx::shared_future<vector_data> middle_uv_data =
        middle_uv.get_data(CENTER);
    
   /* //do local computation first
     hpx::future<scalar_data> next_temperature_middle = 
        hpx::dataflow(
            hpx::launch::async,
            hpx::util::unwrapped(
                [middle_temperature, flag_data, re, pr, dx, dy, dt, alpha]
                (scalar_data const& m_temp, vector_data const& m_uv)
                -> scalar_data
                {
                    uint size_x = m_temp.size_x();
                    uint size_y = m_temp.size_y();

                    scalar_data next(size_x, size_y);

                    for (uint j = 1; j < size_y - 1; j++)
                        for (uint i = 1; i < size_x - 1; i++)
                            compute_temperature_for_cell(
                                    next.get_cell_ref(i, j),
                                    m_temp.get_cell(i, j),
                                    m_temp.get_cell(i - 1, j),
                                    m_temp.get_cell(i + 1, j),
                                    m_temp.get_cell(i, j - 1),
                                    m_temp.get_cell(i, j + 1),
                                    m_uv.get_cell(i, j),
                                    m_uv.get_cell(i - 1, j),
                                    m_uv.get_cell(i, j - 1),
                                    flag_data[j * size_x + i],
                                    re, pr, dx, dy, dt, alpha);

                    return next;
                }
            ),
            middle_temperature_data,
            middle_uv_data
        );*/

            
    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [middle_temperature, flag_data, re, pr, dx, dy, dt, alpha]
            (/*scalar_data next,*/ scalar_data const& m_temp, scalar_data const& l_temp,
            scalar_data const& r_temp, scalar_data const& b_temp, scalar_data const& t_temp,
            vector_data const& m_uv, vector_data const& l_uv,
            vector_data const& b_uv)
            -> scalar_partition
            {
                uint size_x = m_temp.size_x();
                uint size_y = m_temp.size_y();
                   
                scalar_data next(m_temp);

                    for (uint j = 0; j < size_y; j++)
                        for (uint i = 0; i < size_x; i++)
                            compute_temperature_for_cell(
                                next.get_cell_ref(i, j),
                                m_temp.get_cell(i, j),
                                get_left_neighbor(m_temp, l_temp, i, j),
                                get_right_neighbor(m_temp, r_temp, i, j),
                                get_bottom_neighbor(m_temp, b_temp, i, j),
                                get_top_neighbor(m_temp, t_temp, i, j),
                                m_uv.get_cell(i, j),
                                get_left_neighbor(m_uv, l_uv, i, j),
                                get_bottom_neighbor(m_uv, b_uv, i, j),
                                flag_data[j * size_x + i],
                                re, pr, dx, dy, dt, alpha);              
               /* //left and right
                for (uint j = 0; j < size_y; j++)
                {
                    uint i = 0;                     
                    compute_temperature_for_cell(
                        next.get_cell_ref(i, j),
                        m_temp.get_cell(i, j),
                        get_left_neighbor(m_temp, l_temp, i, j),
                        get_right_neighbor(m_temp, r_temp, i, j),
                        get_bottom_neighbor(m_temp, b_temp, i, j),
                        get_top_neighbor(m_temp, t_temp, i, j),
                        m_uv.get_cell(i, j),
                        get_left_neighbor(m_uv, l_uv, i, j),
                        get_bottom_neighbor(m_uv, b_uv, i, j),
                        flag_data[j * size_x + i],
                        re, pr, dx, dy, dt, alpha);
                     
                    
                    i = size_x - 1;
                    compute_temperature_for_cell(
                        next.get_cell_ref(i, j),
                        m_temp.get_cell(i, j),
                        get_left_neighbor(m_temp, l_temp, i, j),
                        get_right_neighbor(m_temp, r_temp, i, j),
                        get_bottom_neighbor(m_temp, b_temp, i, j),
                        get_top_neighbor(m_temp, t_temp, i, j),
                        m_uv.get_cell(i, j),
                        get_left_neighbor(m_uv, l_uv, i, j),
                        get_bottom_neighbor(m_uv, b_uv, i, j),
                        flag_data[j * size_x + i],
                        re, pr, dx, dy, dt, alpha);
                }    
                
                //bottom and top
                for (uint i = 0; i < size_x; i++)
                {
                    uint j = 0;
                    compute_temperature_for_cell(
                        next.get_cell_ref(i, j),
                        m_temp.get_cell(i, j),
                        get_left_neighbor(m_temp, l_temp, i, j),
                        get_right_neighbor(m_temp, r_temp, i, j),
                        get_bottom_neighbor(m_temp, b_temp, i, j),
                        get_top_neighbor(m_temp, t_temp, i, j),
                        m_uv.get_cell(i, j),
                        get_left_neighbor(m_uv, l_uv, i, j),
                        get_bottom_neighbor(m_uv, b_uv, i, j),
                        flag_data[j * size_x + i],
                        re, pr, dx, dy, dt, alpha);
                     
                    
                    j = size_y - 1;
                    compute_temperature_for_cell(
                        next.get_cell_ref(i, j),
                        m_temp.get_cell(i, j),
                        get_left_neighbor(m_temp, l_temp, i, j),
                        get_right_neighbor(m_temp, r_temp, i, j),
                        get_bottom_neighbor(m_temp, b_temp, i, j),
                        get_top_neighbor(m_temp, t_temp, i, j),
                        m_uv.get_cell(i, j),
                        get_left_neighbor(m_uv, l_uv, i, j),
                        get_bottom_neighbor(m_uv, b_uv, i, j),
                        flag_data[j * size_x + i],
                        re, pr, dx, dy, dt, alpha);
                }                     */
                          

                return scalar_partition(middle_temperature.get_id(), next);
            }
        ),
        //std::move(next_temperature_middle),
        middle_temperature_data,
        left_temperature.get_data(LEFT),
        right_temperature.get_data(RIGHT),
        bottom_temperature.get_data(BOTTOM),
        top_temperature.get_data(TOP),
        middle_uv_data,
        left_uv.get_data(LEFT),
        bottom_uv.get_data(BOTTOM)
    );        
}

void custom_grain_size::compute_temperature_for_cell(scalar_cell& middle_temperature,
    scalar_cell const& old_middle_temperature,
    scalar_cell const& left_temperature, scalar_cell const& right_temperature,
    scalar_cell const& bottom_temperature, scalar_cell const& top_temperature,
    vector_cell const& middle_uv, vector_cell const& left_uv,
    vector_cell const& bottom_uv, std::bitset<5> const& type,
    RealType re, RealType pr, RealType dx, RealType dy, RealType dt,
    RealType alpha
    )
{
    
    if (type.test(4))
    {        
        middle_temperature.value =
            dt * (1./re*1./pr * (second_derivative_fwd_bkwd_x(right_temperature.value, old_middle_temperature.value, left_temperature.value, dx)
                                                + second_derivative_fwd_bkwd_y(top_temperature.value, old_middle_temperature.value, bottom_temperature.value, dy)
                                                )
                               // + (i == 1 ? 1 : 0)
                                - first_derivative_u_temp_x(middle_uv.first, left_uv.first, old_middle_temperature.value, right_temperature.value, left_temperature.value, dx, alpha)
                                - first_derivative_v_temp_y(middle_uv.second, bottom_uv.second, old_middle_temperature.value, top_temperature.value, bottom_temperature.value, dy, alpha)
                                )
                                + old_middle_temperature.value;
    }
    else
    {
        middle_temperature.value = old_middle_temperature.value;
    }
}

scalar_partition custom_grain_size::compute_right_hand_side_on_fluid_cells(
    vector_partition const& middle_fg,  vector_partition const& left_fg,
    vector_partition const& bottom_fg,
    std::vector<std::bitset<5> > const& flag_data,
    RealType dx, RealType dy, RealType dt)
{
    hpx::shared_future<vector_data> middle_fg_data = middle_fg.get_data(CENTER);
    
    //do local computation first
    hpx::future<scalar_data> next_rhs_middle = 
        hpx::dataflow(
            hpx::launch::async,
            hpx::util::unwrapped(
                [middle_fg, flag_data, dx, dy, dt]
                (vector_data const& m_fg)
                -> scalar_data
                {
                    uint size_x = m_fg.size_x();
                    uint size_y = m_fg.size_y();

                    scalar_data next(size_x, size_y);

                    for (uint j = 1; j < size_y - 1; j++)
                        for (uint i = 1; i < size_x - 1; i++)
                            compute_rhs_for_cell(
                                next.get_cell_ref(i, j),
                                m_fg.get_cell(i, j),
                                m_fg.get_cell(i - 1, j),
                                m_fg.get_cell(i, j - 1),
                                flag_data[j * size_x + i],
                                dx, dy, dt);

                    return next;
                }
            ),
            middle_fg_data
        );

    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [middle_fg, flag_data, dx, dy, dt]
            (scalar_data next, vector_data const& m_fg, vector_data const& l_fg,
                vector_data const& b_fg)
            -> scalar_partition
            {
                uint size_x = next.size_x();
                uint size_y = next.size_y();
                                
                //left and right
                for (uint j = 0; j < size_y; j++)
                {
                    uint i = 0;                     
                    compute_rhs_for_cell(
                        next.get_cell_ref(i, j),
                        m_fg.get_cell(i, j),
                        get_left_neighbor(m_fg, l_fg, i, j),
                        get_bottom_neighbor(m_fg, b_fg, i, j),
                        flag_data[j * size_x + i],
                        dx, dy, dt);
                     
                    
                    i = size_x - 1;
                    compute_rhs_for_cell(
                        next.get_cell_ref(i, j),
                        m_fg.get_cell(i, j),
                        get_left_neighbor(m_fg, l_fg, i, j),
                        get_bottom_neighbor(m_fg, b_fg, i, j),
                        flag_data[j * size_x + i],
                        dx, dy, dt);
                }    
                
                //bottom and top
                for (uint i = 0; i < size_x; i++)
                {
                    uint j = 0;
                    compute_rhs_for_cell(
                        next.get_cell_ref(i, j),
                        m_fg.get_cell(i, j),
                        get_left_neighbor(m_fg, l_fg, i, j),
                        get_bottom_neighbor(m_fg, b_fg, i, j),
                        flag_data[j * size_x + i],
                        dx, dy, dt);
                     
                    
                    j = size_y - 1;
                    compute_rhs_for_cell(
                        next.get_cell_ref(i, j),
                        m_fg.get_cell(i, j),
                        get_left_neighbor(m_fg, l_fg, i, j),
                        get_bottom_neighbor(m_fg, b_fg, i, j),
                        flag_data[j * size_x + i],
                        dx, dy, dt);
                }                     
                
                return scalar_partition(middle_fg.get_id(), next);
            }
        ),
        std::move(next_rhs_middle),
        middle_fg_data,
        left_fg.get_data(LEFT),
        bottom_fg.get_data(BOTTOM)
    );    
}

void custom_grain_size::compute_rhs_for_cell(scalar_cell& middle_rhs,
    vector_cell const& middle_fg, vector_cell const& left_fg,
    vector_cell const& bottom_fg, std::bitset<5> const& type,
    RealType dx, RealType dy, RealType dt)
{
    if (type.test(4))          
        middle_rhs.value =
            1./dt * ( (middle_fg.first - left_fg.first)/dx
                     + (middle_fg.second - bottom_fg.second)/dy);
}

scalar_partition custom_grain_size::set_pressure_on_boundary_and_obstacles(
        scalar_partition const& middle_p, scalar_partition const& left_p,
        scalar_partition const& right_p, scalar_partition const& bottom_p,
        scalar_partition const& top_p, std::vector<std::bitset<5> > const& flag_data)
{ 
    //do local computation first
    hpx::future<scalar_data> next_p_middle = 
        hpx::dataflow(
            hpx::launch::async,
            hpx::util::unwrapped(
                [middle_p, flag_data]
                (scalar_data m_p)
                -> scalar_data
                {
                    uint size_x = m_p.size_x();
                    uint size_y = m_p.size_y();

                    //scalar_data next(size_x, size_y);

                    for (uint j = 1; j < size_y - 1; j++)
                        for (uint i = 1; i < size_x - 1; i++)
                            set_pressure_for_cell(
                                m_p.get_cell_ref(i, j),
                                m_p.get_cell(i - 1, j),
                                m_p.get_cell(i + 1, j),
                                m_p.get_cell(i, j - 1),
                                m_p.get_cell(i, j + 1),
                                flag_data[j * size_x + i]);

                    return m_p;
                }
            ),
            middle_p.get_data(CENTER)
        );
   
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
                             
                //left and right
                for (uint j = 0; j < size_y; j++)
                {
                    uint i = 0;                     
                    set_pressure_for_cell(
                       next.get_cell_ref(i, j),
                       get_left_neighbor(next, l_p, i, j),
                       get_right_neighbor(next, r_p, i, j),
                       get_bottom_neighbor(next, b_p, i, j),
                       get_top_neighbor(next, t_p, i, j),
                       flag_data[j * size_x + i]);
                     
                    
                    i = size_x - 1;
                    set_pressure_for_cell(
                       next.get_cell_ref(i, j),
                       get_left_neighbor(next, l_p, i, j),
                       get_right_neighbor(next, r_p, i, j),
                       get_bottom_neighbor(next, b_p, i, j),
                       get_top_neighbor(next, t_p, i, j),
                       flag_data[j * size_x + i]);
                }    
                
                //bottom and top
                for (uint i = 1; i < size_x - 1; i++)
                {
                    uint j = 0;
                    set_pressure_for_cell(
                       next.get_cell_ref(i, j),
                       get_left_neighbor(next, l_p, i, j),
                       get_right_neighbor(next, r_p, i, j),
                       get_bottom_neighbor(next, b_p, i, j),
                       get_top_neighbor(next, t_p, i, j),
                       flag_data[j * size_x + i]);
                     
                    
                    j = size_y - 1;
                    set_pressure_for_cell(
                       next.get_cell_ref(i, j),
                       get_left_neighbor(next, l_p, i, j),
                       get_right_neighbor(next, r_p, i, j),
                       get_bottom_neighbor(next, b_p, i, j),
                       get_top_neighbor(next, t_p, i, j),
                       flag_data[j * size_x + i]);
                }                                              
                
                return scalar_partition(middle_p.get_id(), next);
            }
        ),
        std::move(next_p_middle),
        left_p.get_data(LEFT),
        right_p.get_data(RIGHT),
        bottom_p.get_data(BOTTOM),
        top_p.get_data(TOP)
    );
}

void custom_grain_size::set_pressure_for_cell(scalar_cell& middle_p,
    scalar_cell const& left_p, scalar_cell const& right_p,
    scalar_cell const& bottom_p, scalar_cell const& top_p,
    std::bitset<5> const& type)
{
    //east
    if (type == std::bitset<5>("01000") || type == std::bitset<5>("00111"))
        middle_p.value = right_p.value;

    //west
    else if (type == std::bitset<5>("00100") || type == std::bitset<5>("01011"))
        middle_p.value = left_p.value;

    //south
    else if (type == std::bitset<5>("00010") || type == std::bitset<5>("01101"))
        middle_p.value = bottom_p.value;

    //north
    else if (type == std::bitset<5>("00001") || type == std::bitset<5>("01110"))
        middle_p.value = top_p.value;

    //NE
    else if (type == std::bitset<5>("01001"))
        middle_p.value = (top_p.value + right_p.value) / 2.;

    //SE
    else if (type == std::bitset<5>("01010"))
        middle_p.value = (bottom_p.value + right_p.value) / 2.;

    //SW
    else if (type == std::bitset<5>("00110"))
        middle_p.value = (bottom_p.value + left_p.value) / 2.;

    //NW
    else if (type == std::bitset<5>("00101"))
        middle_p.value = (top_p.value + left_p.value) / 2.; 
     
    
}

scalar_partition custom_grain_size::sor_cycle(scalar_partition const& middle_p,
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

                //scalar_data next(size_x, size_y);
                for (uint j = 0; j < size_y - 1; j++)
                    for (uint i = 0; i < size_x - 1; i++)
                    {
                        auto& top_type =
                            flag_data[(j + 1) * size_x + i];
                        auto& right_type =
                            flag_data[j * size_x + i + 1];
                     
                        auto& type = flag_data[j * size_x + i];

                        //top is obstacle or boundary cell
                        if (!top_type.test(4))
                            set_pressure_for_cell(
                                m_p.get_cell_ref(i, j + 1),
                                get_left_neighbor(m_p, l_p, i, j + 1),
                                m_p.get_cell(i + 1, j + 1),
                                m_p.get_cell(i, j),
                                get_top_neighbor(m_p, t_p, i, j + 1),
                                top_type);

                        if (!right_type.test(4))
                            set_pressure_for_cell(
                                m_p.get_cell_ref(i + 1, j),
                                m_p.get_cell(i, j),
                                get_right_neighbor(m_p, r_p, i + 1, j),
                                get_bottom_neighbor(m_p, b_p, i + 1, j),
                                m_p.get_cell(i + 1, j + 1),
                                right_type);

                        //current cell is fluid cell

                        if (type.test(4))
                            do_sor_cycle_for_cell(
                                m_p.get_cell_ref(i, j),
                                get_left_neighbor(m_p, l_p, i, j),
                                m_p.get_cell(i + 1, j),
                                get_bottom_neighbor(m_p, b_p, i, j),
                                m_p.get_cell(i, j + 1),
                                m_rhs.get_cell(i, j),
                                flag_data[j * size_x + i],
                                dx_sq, dy_sq, part1, part2);
                        else if (i == 0 && j == 0)
                            set_pressure_for_cell(
                                m_p.get_cell_ref(i, j),
                                get_left_neighbor(m_p, l_p, i, j),
                                m_p.get_cell(i + 1, j),
                                get_bottom_neighbor(m_p, b_p, i, j),
                                m_p.get_cell(i, j + 1),
                                flag_data[j * size_x + i]);
                    }
                
                //TODO: get proper value from right and top for remote partition
                //right
                uint i = size_x - 1;
                for (uint j = 0; j < size_y - 1; j++)
                {
                    do_sor_cycle_for_cell(
                        m_p.get_cell_ref(i, j),
                        m_p.get_cell(i - 1, j),
                        get_right_neighbor(m_p, r_p, i, j),
                        get_bottom_neighbor(m_p, b_p, i, j),
                        m_p.get_cell(i, j + 1),
                        m_rhs.get_cell(i, j),
                        flag_data[j * size_x + i],
                        dx_sq, dy_sq, part1, part2);
                }
                
                //top
                uint j = size_y - 1;
                for (i = 0; i < size_x; i++)
                {                                       
                    do_sor_cycle_for_cell(
                        m_p.get_cell_ref(i, j),
                        get_left_neighbor(m_p, l_p, i, j),
                        get_right_neighbor(m_p, r_p, i, j),
                        m_p.get_cell(i, j - 1),
                        get_top_neighbor(m_p, t_p, i, j),
                        m_rhs.get_cell(i, j),
                        flag_data[j * size_x + i],
                        dx_sq, dy_sq, part1, part2);                  
                    }                
               
                i = size_x - 1;
                j = size_y - 1;
                
                set_pressure_for_cell(
                                m_p.get_cell_ref(i, j),
                                m_p.get_cell(i - 1, j),
                                get_right_neighbor(m_p, r_p, i, j),
                                m_p.get_cell(i, j - 1),
                                get_top_neighbor(m_p, t_p, i, j),
                                flag_data[j * size_x + i]);
                
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

void custom_grain_size::do_sor_cycle_for_cell(scalar_cell& middle_p,
    scalar_cell const& left_p, scalar_cell const& right_p,
    scalar_cell const& bottom_p, scalar_cell const& top_p,
    scalar_cell const& middle_rhs, std::bitset<5> const& type,
    RealType dx_sq, RealType dy_sq, RealType part1, RealType part2)
{    
    if (type.test(4))
                middle_p.value = part1 * middle_p.value
                            + part2 * ( (right_p.value + left_p.value) / dx_sq 
                                        + (top_p.value + bottom_p.value) / dy_sq 
                                        - middle_rhs.value);
}

hpx::future<RealType> custom_grain_size::compute_residual(scalar_partition const& middle_p,
                            scalar_partition const& left_p,
                            scalar_partition const& right_p,
                            scalar_partition const& bottom_p,
                            scalar_partition const& top_p,
                            scalar_partition const& middle_rhs,
                            std::vector<std::bitset<5> > const& flag_data,
                            RealType dx, RealType dy)
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
                
                RealType local_residual = 0;
                
                for (uint j = 1; j < size_y - 1; j++)
                        for (uint i = 1; i < size_x - 1; i++)
                            local_residual +=
                                compute_residual_for_cell(
                                    m_p.get_cell(i, j),
                                    get_left_neighbor(m_p, l_p, i, j),
                                    get_right_neighbor(m_p, r_p, i, j),
                                    get_bottom_neighbor(m_p, b_p, i, j),
                                    get_top_neighbor(m_p, t_p, i, j),
                                    m_rhs.get_cell(i, j),
                                    flag_data[j * size_x + i],
                                    over_dx_sq, over_dy_sq);
                
                return local_residual;
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

RealType custom_grain_size::compute_residual_for_cell(scalar_cell const& middle_p,
    scalar_cell const& left_p, scalar_cell const& right_p,
    scalar_cell const& bottom_p, scalar_cell const& top_p,
    scalar_cell const& middle_rhs, std::bitset<5> const& type,
    RealType over_dx_sq, RealType over_dy_sq)
{
    if (type.test(4))
    {
        RealType tmp =
            (right_p.value - 2*middle_p.value + left_p.value)*over_dx_sq 
            + (top_p.value - 2*middle_p.value + bottom_p.value)*over_dy_sq 
            - middle_rhs.value;
        
        return std::pow(tmp, 2);
    }
    
    return 0;
}

hpx::future<std::pair<vector_partition, std::pair<RealType, RealType> > > 
custom_grain_size::update_velocities(
       vector_partition const& middle_uv, scalar_partition const& middle_p,
       scalar_partition const& right_p, scalar_partition const& top_p, 
       vector_partition const& middle_fg, std::vector<std::bitset<5> > const& flag_data,
       RealType dx, RealType dy, RealType dt)
{
    RealType const over_dx = 1./dx;
    RealType const over_dy = 1./dy;
       
    hpx::shared_future<scalar_data> middle_p_data =
        middle_p.get_data(CENTER);
    
    hpx::shared_future<vector_data> middle_fg_data =
        middle_fg.get_data(CENTER);
    
    //do local computation first
    hpx::future<std::pair<vector_data, std::pair<RealType, RealType> > > next_middle_uv = 
        hpx::dataflow(
            hpx::launch::async,
            hpx::util::unwrapped(
                [middle_uv, flag_data, over_dx, over_dy, dt]
                (vector_data m_uv, scalar_data const& m_p,
                    vector_data const& m_fg)
                -> std::pair<vector_data, std::pair<RealType, RealType> >
                {
                    uint size_x = m_uv.size_x();
                    uint size_y = m_uv.size_y();

                    std::pair<RealType, RealType> max_uv(0., 0.);
                    
                    //scalar_data next(size_x, size_y);
    
                    for (uint j = 0; j < size_y - 1; j++)
                        for (uint i = 0; i < size_x - 1; i++)
                        {
                            vector_cell& middle_cell = m_uv.get_cell_ref(i, j);
                            
                            auto& type = flag_data[j * size_x + i];
                            
                            update_velocity_for_cell(
                                middle_cell,
                                m_p.get_cell(i, j),
                                m_p.get_cell(i + 1, j),
                                m_p.get_cell(i, j + 1),
                                m_fg.get_cell(i, j),
                                type,
                                over_dx, over_dy, dt);
                            
                            if (type.test(4))
                            {
                                max_uv.first = (std::abs(middle_cell.first) > max_uv.first)
                                                ? std::abs(middle_cell.first) : max_uv.first;

                                max_uv.second = (std::abs(middle_cell.second) > max_uv.second)
                                                ? std::abs(middle_cell.second) : max_uv.second; 
                            }
                        }

                    return std::make_pair(m_uv, max_uv);
                }
            ),
            middle_uv.get_data(CENTER),
            middle_p_data,
            middle_fg_data
        );
    
    return hpx::dataflow(
        hpx::launch::async,
        hpx::util::unwrapped(
            [middle_uv, flag_data, over_dx, over_dy, dt]
            (std::pair<vector_data, std::pair<RealType, RealType> > next_, scalar_data const& m_p, scalar_data const& r_p,
                scalar_data const& t_p, vector_data const& m_fg)
            -> std::pair<vector_partition, std::pair<RealType, RealType> >
            {
                vector_data next = next_.first;
                auto max_uv = next_.second;
                
                uint size_x = next.size_x();
                uint size_y = next.size_y();

                //right
                for (uint j = 0; j < size_y; j++)
                {
                    uint i = size_x - 1;
                    auto& type = flag_data[j * size_x + i];
                    
                    vector_cell& middle_cell = next.get_cell_ref(i, j);
                    
                    update_velocity_for_cell(
                        next.get_cell_ref(i, j),
                        m_p.get_cell(i, j),
                        get_right_neighbor(m_p, r_p, i, j),
                        get_top_neighbor(m_p, t_p, i, j),
                        m_fg.get_cell(i, j),
                        type,
                        over_dx, over_dy, dt);
                    
                    if (type.test(4))
                    {
                        max_uv.first = (std::abs(middle_cell.first) > max_uv.first)
                                        ? std::abs(middle_cell.first) : max_uv.first;

                        max_uv.second = (std::abs(middle_cell.second) > max_uv.second)
                                        ? std::abs(middle_cell.second) : max_uv.second; 
                    }        
                }
                
                //top
                for (uint i = 0; i < size_x -1; i++)
                {
                    uint j = size_y - 1;
                    auto& type = flag_data[j * size_x + i];
 
                    vector_cell& middle_cell = next.get_cell_ref(i, j);

                    update_velocity_for_cell(
                        next.get_cell_ref(i, j),
                        m_p.get_cell(i, j),
                        get_right_neighbor(m_p, r_p, i, j),
                        get_top_neighbor(m_p, t_p, i, j),
                        m_fg.get_cell(i, j),
                        type,
                        over_dx, over_dy, dt);
                    
                    if (type.test(4))
                           {
                               max_uv.first = (std::abs(middle_cell.first) > max_uv.first)
                                               ? std::abs(middle_cell.first) : max_uv.first;

                               max_uv.second = (std::abs(middle_cell.second) > max_uv.second)
                                               ? std::abs(middle_cell.second) : max_uv.second; 
                           }
                }
                
                return std::make_pair(vector_partition(middle_uv.get_id(), next), max_uv);
            }
        ),
        next_middle_uv,
        middle_p_data,
        right_p.get_data(RIGHT),
        top_p.get_data(TOP),
        middle_fg_data
    );    
}

void custom_grain_size::update_velocity_for_cell(vector_cell& middle_uv,
        scalar_cell const& middle_p, scalar_cell const& right_p,
        scalar_cell const& top_p, vector_cell const& middle_fg,
        std::bitset<5> const& type, RealType over_dx, RealType over_dy, RealType dt)
{
    if (type.test(4))
    {
        if (type.test(3))
            middle_uv.first = middle_fg.first - dt * over_dx * (right_p.value - middle_p.value);

        if (type.test(0))
            middle_uv.second = middle_fg.second - dt * over_dy * (top_p.value - middle_p.value);
    }
}

void custom_grain_size::compute_stream_vorticity_heat(scalar_data& stream_center, scalar_data& vorticity_center, scalar_data& heat_center,
                                                    scalar_data const& stream_bottom, scalar_data const& heat_bottom,
                                                    vector_data const& uv_center, vector_data const& uv_right, vector_data const& uv_top,
                                                    scalar_data const& temp_center, scalar_data const& temp_right,
                                                    std::vector<std::bitset<5> > const& flag_data,
                                                    uint global_i, uint global_j, uint i_max, uint j_max, RealType re, RealType pr,
                                                    RealType dx, RealType dy)
{
    uint size_x = stream_center.size_x();
    uint size_y = stream_center.size_y();

    for (uint j = 0; j < size_y; j++)
        for (uint i = 0; i < size_x; i++)
        {
            std::bitset<5> cell_type = flag_data[j*size_x + i];


            if (in_range(0, i_max, 1, j_max, global_i + i, global_j + j))
            {
                if (cell_type.test(4) || global_i + i == 0)
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

