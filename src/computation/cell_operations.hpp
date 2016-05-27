/** Header that declares all the necessary methods for computation on the grids
 */

#ifndef CELL_OPERATIONS_HPP
#define CELL_OPERATIONS_HPP

#include <hpx/hpx.hpp>

#include "grid/types.hpp"
#include "util/typedefs.hpp"
#include "util/boundary_data.hpp"
#include "stencils.hpp"

namespace computation
{ 
 
inline void set_velocity_for_obstacle(vector_cell& middle,
    vector_cell  left, vector_cell  right,
    vector_cell  bottom, vector_cell  top,
    std::bitset<5>  cell_type)
{
    middle.first = -bottom.first * cell_type.test(has_fluid_south)
                    -top.first * cell_type.test(has_fluid_north);
                    
    middle.second = -left.second * cell_type.test(has_fluid_west)
                    -right.second * cell_type.test(has_fluid_east);    
}

inline void set_velocity_for_left_boundary(vector_cell& middle,
    vector_cell  right, int type, Real u_left, Real v_left)
{
    switch(type)
    {
    case noslip:
            middle.first = 0;
            middle.second = -right.second;
            break;
    case slip:
            middle.first = 0;
            middle.second = right.second;
            break;
    case outstream:
            middle.first = right.first;
            middle.second = right.second;
            break;
    case instream:
            middle.first = u_left;
            middle.second = 2 * v_left - right.second;
            break;  
    }
}

inline void set_velocity_for_right_boundary(vector_cell& middle,
    vector_cell& left, vector_cell  left_left,
    int type, Real u_right, Real v_right)
{
    switch(type)
    {
    case noslip:
            left.first = 0;
            middle.second = -left.second;
            break;
    case slip:
            left.first = 0;
            middle.second = left.second;
            break;
    case outstream:
            left.first = left_left.first;
            middle.second = left.second;                
            break;
    case instream:
            left.first = u_right;
            middle.second = 2 * v_right - left.second;
            break;  
    }      
}

inline void set_velocity_for_bottom_boundary(vector_cell& middle,
    vector_cell  top, int type, Real u_bottom, Real v_bottom)
{
    switch(type)
    {
    case noslip:
            middle.first = -top.first;
            middle.second = 0;
            break;
    case slip:
            middle.first = top.first;
            middle.second = 0;
            break;
    case outstream:
            middle.first = top.first;
            middle.second = top.second;
            break;
    case instream:
            middle.first = 2 * u_bottom - top.first;
            middle.second = v_bottom;
            break;  
    }
}

inline void set_velocity_for_top_boundary(vector_cell& middle,
    vector_cell& bottom, vector_cell  bottom_bottom,
    int type, Real u_top, Real v_top)
{
    switch(type)
    {
    case noslip:
            middle.first = 2 * u_top - bottom.first;
            bottom.second = 0;
            break;
    case slip:
            middle.first = bottom.first;
            bottom.second = 0;
            break;
    case outstream:
            middle.first = bottom.first;
            bottom.second = bottom_bottom.second;
            break;
    case instream:
            middle.first = 2 * u_top - bottom.first;
            bottom.second = v_top;
            break;  
    }  
}
    
inline void set_velocity_for_fluid(vector_cell& middle,
    vector_cell  left, vector_cell  right, 
    vector_cell  bottom, vector_cell  top,
    std::bitset<5>  cell_type)
{
    if (!cell_type.test(has_fluid_north))
        middle.second = 0;
        
    if (!cell_type.test(has_fluid_east))
        middle.first = 0;
}   

inline void set_temperature_for_boundary(Real& middle,
    Real  neighbor, int type, Real temperature_boundary, uint i,
    Real dx, Real dy)
{
    switch(type)
    {
    case noslip: 
            middle = 2 * temperature_boundary - neighbor ;
            break;            
    case slip:
            middle = neighbor + dx * temperature_boundary* ((i - 0.5) * dy);
            break;
    }
}

inline void compute_fg_for_cell(vector_cell& middle_fg,
        vector_cell  middle_uv, vector_cell  left_uv,
        vector_cell  right_uv, vector_cell  bottom_uv,
        vector_cell  top_uv, vector_cell  bottomright_uv,
        vector_cell  topleft_uv, Real  middle_temperature,
        Real  right_temperature,
        Real  top_temperature, std::bitset<5>  type,
        Real re, Real gx, Real gy, Real beta,
        Real dx, Real dy, Real dt, Real alpha)
{      
    middle_fg.first =
        middle_uv.first
        + type.test(3) * (dt * (
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
                              * (middle_temperature 
                                    + right_temperature )
                              * gx
                        );
    middle_fg.second =
        middle_uv.second
            + type.test(0) * (dt * (
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
                                * (middle_temperature 
                                    + top_temperature )
                                * gy
                            );
}

inline Real compute_temperature_for_cell(Real middle_temperature,
    Real left_temperature, Real right_temperature,
    Real bottom_temperature, Real top_temperature,
    vector_cell  middle_uv, vector_cell  left_uv,
    vector_cell  bottom_uv, Real re, Real pr, Real dx,
    Real dy, Real dt, Real alpha)
{         
   return dt * (
            1./re*1./pr * (
                    second_derivative_fwd_bkwd_x(right_temperature ,
                        middle_temperature ,
                        left_temperature , dx)
        
                    + second_derivative_fwd_bkwd_y(top_temperature ,
                        middle_temperature ,
                        bottom_temperature  ,dy)
            )
                
                           // + (i == 1 ? 1 : 0)
            - first_derivative_u_temp_x(middle_uv.first, left_uv.first,
                middle_temperature , right_temperature ,
                left_temperature , dx, alpha)

            - first_derivative_v_temp_y(middle_uv.second,
                bottom_uv.second, middle_temperature ,
                top_temperature , bottom_temperature , dy,
                alpha)
            )
            + middle_temperature ;
}

inline Real compute_rhs_for_cell(
    vector_cell  middle_fg, vector_cell  left_fg,
    vector_cell  bottom_fg, Real dx, Real dy, Real dt)
{        
    return 1./dt * ( (middle_fg.first - left_fg.first)/dx
                     + (middle_fg.second - bottom_fg.second)/dy);
}


inline void set_pressure_for_cell(Real& middle_p,
    Real left_p, Real right_p, Real bottom_p, Real top_p,
    std::bitset<5>  type)
    {
        middle_p  = (left_p * type.test(has_fluid_west)
                        + right_p * type.test(has_fluid_east)
                        + bottom_p * type.test(has_fluid_south)
                        + top_p * type.test(has_fluid_north))
                        / (type.test(has_fluid_west) + type.test(has_fluid_east)
                            + type.test(has_fluid_south) 
                            + type.test(has_fluid_north)
                        );        
    }
    
inline void set_pressure_for_boundary(Real& middle_p,
    Real left_p, Real right_p, Real bottom_p, Real top_p,
    std::bitset<5>  type)
{
    middle_p  = left_p * !type.test(has_fluid_west)
                    + right_p * !type.test(has_fluid_east)
                    + bottom_p * !type.test(has_fluid_south)
                    + top_p * !type.test(has_fluid_north);        
}
    
inline Real do_sor_cycle_for_cell(Real middle_p,
    Real left_p, Real right_p, Real bottom_p, Real top_p,
    Real middle_rhs, Real dx_sq, Real dy_sq, Real part1,
    Real part2)
{   
    return part1 * middle_p 
            + part2 * (
                    (right_p + left_p ) / dx_sq 
                    + (top_p + bottom_p ) / dy_sq 
                    - middle_rhs 
                  );
}

inline Real compute_residual_for_cell(Real middle_p,
    Real left_p, Real right_p, Real bottom_p, Real top_p,
    Real middle_rhs, Real over_dx_sq, Real over_dy_sq)
{
    Real tmp =
        (right_p - 2 * middle_p + left_p ) * over_dx_sq 
        + (top_p - 2 * middle_p + bottom_p ) * over_dy_sq 
        - middle_rhs;
    
    return std::pow(tmp, 2);
}



inline void update_velocity_for_cell(vector_cell& middle_uv,
    Real middle_p, Real right_p, Real top_p,
    vector_cell  middle_fg, std::bitset<5>  type, Real over_dx,
    Real over_dy, Real dt)
{
    if (type.test(3))
        middle_uv.first =
            middle_fg.first - dt * over_dx
            * (right_p  - middle_p );

    if (type.test(0))
        middle_uv.second =
            middle_fg.second - dt * over_dy
            * (top_p  - middle_p );
}
}//end namespace computation

#endif /* CELL_OPERATIONS_HPP */

