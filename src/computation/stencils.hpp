/** All stencils are provided in this header as inline methods
 */
#ifndef COMPUTATION_STENCILS_HPP
#define COMPUTATION_STENCILS_HPP

#include "util/typedefs.hpp"
#include <cmath>
#include <cstdlib>

namespace computation {

inline RealType second_derivative_fwd_bkwd_x(
    RealType right, RealType middle, RealType left, RealType dx)
{
    return (right - 2 * middle + left) / std::pow(dx, 2);
}

inline RealType second_derivative_fwd_bkwd_y(
    RealType top, RealType middle, RealType bottom, RealType dy)
{
    return (top - 2 * middle + bottom) / std::pow(dy, 2);
}

inline RealType first_derivative_fwd_x(
    RealType right, RealType middle, RealType dx)
{
    return (right - middle) / dx;
}

inline RealType first_derivative_fwd_y(
    RealType top, RealType middle, RealType dy)
{
    return (top - middle) / dy;
}

inline RealType first_derivative_bkwd_x(
    RealType middle, RealType left, RealType dx)
{
    return (middle - left) / dx;
}

inline RealType first_derivative_bkwd_y(
    RealType middle, RealType bottom, RealType dy)
{
    return (middle - bottom) / dy;
}

inline RealType first_derivative_of_square_x(
    RealType right, RealType middle, RealType left, RealType dx,
    RealType alpha = 0.9)
{
    return 1./dx * (std::pow((middle + right) / 2., 2)
                    - std::pow((left + middle) / 2., 2))
            + alpha / dx
            * (std::abs(middle + right) * (middle - right) / 4.
                - std::abs(left + middle) * (left - middle) / 4.);
}

inline RealType first_derivative_of_square_y(
    RealType top, RealType middle, RealType bottom, RealType dy,
    RealType alpha = 0.9)
{
    return 1./dy * (std::pow((middle + top) / 2., 2)
                        - std::pow((bottom + middle) / 2., 2))
            + alpha / dy
            * (std::abs(middle + top) * (middle - top) / 4.
                - std::abs(bottom + middle) * (bottom - middle) / 4.);
}

inline RealType first_derivative_of_product_x(
    RealType u_left, RealType u_middle, RealType u_top, RealType u_topleft,
    RealType v_left, RealType v_middle, RealType v_right, RealType dx,
    RealType alpha = 0.9)
{
    return 1./dx * ((u_middle + u_top)
                        * (v_middle + v_right) / 4.
                        - (u_left + u_topleft) * (v_left + v_middle) / 4.)
            + alpha / dx * (std::abs(u_middle + u_top)
                * (v_middle - v_right) / 4.
            - std::abs(u_left + u_topleft) * (v_left - v_middle) / 4.);
}

inline RealType first_derivative_of_product_y(
    RealType v_right, RealType v_middle, RealType v_bottom,
    RealType v_bottomright, RealType u_bottom, RealType u_middle,
    RealType u_top, RealType dy, RealType alpha = 0.9)
{
    return 1./dy * ((v_middle + v_right)  * (u_middle + u_top) / 4. 
            - (v_bottom + v_bottomright) * (u_bottom + u_middle) / 4.)
            + alpha / dy * (std::abs(v_middle + v_right) 
                * (u_middle - u_top) / 4.
            - std::abs(v_bottom + v_bottomright) * (u_bottom - u_middle) / 4.);
}

inline RealType first_derivative_u_temp_x(
    RealType u_center, RealType u_left, RealType t_center, RealType t_right,
    RealType t_left, RealType dx, RealType alpha = 0.9)
{
    return 1./dx * (u_center * (t_center + t_right) / 2. 
                - u_left * (t_left + t_center) / 2.)
            + alpha / dx * (std::abs(u_center)*(t_center - t_right) / 2.
            - std::abs(u_left) * (t_left - t_center) / 2.);
}

inline RealType first_derivative_v_temp_y(
    RealType v_center, RealType v_bottom, RealType t_center, RealType t_top, 
    RealType t_bottom, RealType dy, RealType alpha = 0.9)
{
    return 1./dy * (v_center * (t_center + t_top)/ 2.
                - v_bottom * (t_bottom + t_center)/2.)
            + alpha / dy
                * (std::abs(v_center)*(t_center - t_top)/2.
                - std::abs(v_bottom)*(t_bottom - t_center)/2.);
}


}//namespace computation
#endif
