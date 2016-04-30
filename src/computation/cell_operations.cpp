#include "cell_operations.hpp"

#include "stencils.hpp"

namespace computation
{
void set_velocity_for_cell(vector_cell& middle,
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
        switch(static_cast<int>(type.left))
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
                middle.second = 2 * v.left - right.second;
                break;  
        }

    //right
    else if (cell_type == std::bitset<5>("01011"))
        switch(static_cast<int>(type.right))
        {
        case 1: middle.second = -left.second;
                break;
        case 2: middle.second = left.second;
                break;
        case 3: middle.second = left.second;
                break;
        case 4: middle.second = 2 * v.right - left.second;
                break;  
        }      

    //bottom
    else if (cell_type == std::bitset<5>("01110"))
        switch(static_cast<int>(type.bottom))
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
        case 4: middle.first = 2 * u.bottom - top.first;
                middle.second = v.bottom;
                break;  
        }  

    //top
    else if (cell_type == std::bitset<5>("01101"))
        switch(static_cast<int>(type.top))
        {
        case 1: middle.first = 2 * u.top - bottom.first;
                break;
        case 2: middle.first = bottom.first;
                break;
        case 3: middle.first = bottom.first;
                break;
        case 4: middle.first = 2 * u.top - bottom.first;
                break;  
        } 
}  

void set_temperature_for_cell(scalar_cell& middle,
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
        switch(static_cast<int>(boundary_data_type.left))
        {
        case 1: middle.value =
                    2 * temperature_boundary_data.left - right.value;
                break;            
        case 2: middle.value =
                    right.value
                    + dx * temperature_boundary_data.left* ((j - 0.5) * dy);
                break;
        }
    
    //right
    if (cell_type == std::bitset<5>("01011"))
        switch(static_cast<int>(boundary_data_type.right))
        {
        case 1: middle.value = 2 * temperature_boundary_data.right - left.value;
                break;            
        case 2: middle.value =
                    left.value
                    + dx * temperature_boundary_data.right * ((j - 0.5) * dy);
                break;
        }

    //bottom
    if (cell_type == std::bitset<5>("01110"))
        switch(static_cast<int>(boundary_data_type.bottom))
        {
        case 1: middle.value = 2 * temperature_boundary_data.bottom - top.value;
                break;            
        case 2: middle.value =
                    top.value
                    + dy * temperature_boundary_data.bottom * ((i - 0.5) * dx);
                break;
        }

    //top
    if (cell_type == std::bitset<5>("01101"))
        switch(static_cast<int>(boundary_data_type.top))
        {
        case 1: middle.value = 2 * temperature_boundary_data.top - bottom.value;
                break;            
        case 2: middle.value =
                    bottom.value
                + dy * temperature_boundary_data.top * ((i - 0.5) * dx);
                break;
        }    
}

void compute_fg_for_cell(vector_cell& middle_fg,
        vector_cell const& middle_uv, vector_cell const& left_uv,
        vector_cell const& right_uv, vector_cell const& bottom_uv,
        vector_cell const& top_uv, vector_cell const& bottomright_uv,
        vector_cell const& topleft_uv, scalar_cell const& middle_temperature,
        scalar_cell const& right_temperature,
        scalar_cell const& top_temperature, std::bitset<5> const& type,
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
            middle_fg.first =
                middle_uv.first
                + dt * (
                        1./re
                        *   (second_derivative_fwd_bkwd_x(right_uv.first, 
                                middle_uv.first, left_uv.first, dx)
                            + second_derivative_fwd_bkwd_y(top_uv.first,
                                middle_uv.first, bottom_uv.first, dy))

                            - first_derivative_of_square_x(right_uv.first,
                                middle_uv.first, left_uv.first, dx, alpha)
                            - first_derivative_of_product_y(right_uv.second,
                                middle_uv.second, bottom_uv.second,
                                bottomright_uv.second, bottom_uv.first,
                                middle_uv.first, top_uv.first, dy, alpha)

                            + gx
                            )
                            - beta * dt / 2.
                                * (middle_temperature.value
                                    + right_temperature.value)
                                * gx;
        }

        if (!type.test(0))
            middle_fg.second = middle_uv.second;
        else
        {
            middle_fg.second =
                middle_uv.second
                    + dt * (
                            1./re
                            *   (second_derivative_fwd_bkwd_x(right_uv.second,
                                    middle_uv.second, left_uv.second, dx)
                                + second_derivative_fwd_bkwd_y(top_uv.second,
                                    middle_uv.second, bottom_uv.second, dy))

                            - first_derivative_of_product_x(left_uv.first,
                                middle_uv.first, top_uv.first, topleft_uv.first,
                                left_uv.second, middle_uv.second,
                                right_uv.second, dx, alpha)
                
                            - first_derivative_of_square_y(top_uv.second,
                                middle_uv.second, bottom_uv.second, dy, alpha)

                            + gy
                            )
                            - beta * dt / 2.
                                * (middle_temperature.value
                                    + top_temperature.value)
                                * gy;
        }
    }
}


void compute_temperature_for_cell(scalar_cell& middle_temperature,
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
            dt * (
                1./re*1./pr * (
                        second_derivative_fwd_bkwd_x(right_temperature.value,
                            old_middle_temperature.value,
                            left_temperature.value, dx)
            
                        + second_derivative_fwd_bkwd_y(top_temperature.value,
                            old_middle_temperature.value,
                            bottom_temperature.value ,dy)
                )
                    
                               // + (i == 1 ? 1 : 0)
                - first_derivative_u_temp_x(middle_uv.first, left_uv.first,
                    old_middle_temperature.value, right_temperature.value,
                    left_temperature.value, dx, alpha)

                - first_derivative_v_temp_y(middle_uv.second,
                    bottom_uv.second, old_middle_temperature.value,
                    top_temperature.value, bottom_temperature.value, dy,
                    alpha)
                )
                + old_middle_temperature.value;
    }
    else
        middle_temperature.value = old_middle_temperature.value;
}

void compute_rhs_for_cell(scalar_cell& middle_rhs,
    vector_cell const& middle_fg, vector_cell const& left_fg,
    vector_cell const& bottom_fg, std::bitset<5> const& type,
    RealType dx, RealType dy, RealType dt)
{
    if (type.test(4))          
        middle_rhs.value =
            1./dt * ( (middle_fg.first - left_fg.first)/dx
                     + (middle_fg.second - bottom_fg.second)/dy);
}

void set_pressure_for_cell(scalar_cell& middle_p,
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

void do_sor_cycle_for_cell(scalar_cell& middle_p,
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

RealType compute_residual_for_cell(scalar_cell const& middle_p,
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

void update_velocity_for_cell(vector_cell& middle_uv,
    scalar_cell const& middle_p, scalar_cell const& right_p,
    scalar_cell const& top_p, vector_cell const& middle_fg,
    std::bitset<5> const& type, RealType over_dx, RealType over_dy,
    RealType dt)
{
    if (type.test(4))
    {
        if (type.test(3))
            middle_uv.first =
                middle_fg.first - dt * over_dx
                * (right_p.value - middle_p.value);

        if (type.test(0))
            middle_uv.second =
                middle_fg.second - dt * over_dy
                * (top_p.value - middle_p.value);
    }
}

}
