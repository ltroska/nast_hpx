/** All stencils are provided in this header as inline methods
 */
#ifndef NAST_HPX_GRID_FD_STENCILS_HPP
#define NAST_HPX_GRID_FD_STENCILS_HPP

#include "util/typedefs.hpp"
#include <cmath>
#include <cstdlib>

namespace nast_hpx { namespace grid {

inline Real second_derivative_fwd_bkwd_x(
    Real right, Real middle, Real left, Real dx)
{
    return (right - 2 * middle + left) / std::pow(dx, 2);
}

inline Real second_derivative_fwd_bkwd_y(
    Real top, Real middle, Real bottom, Real dy)
{
    return (top - 2 * middle + bottom) / std::pow(dy, 2);
}

inline Real first_derivative_fwd_x(
    Real right, Real middle, Real dx)
{
    return (right - middle) / dx;
}

inline Real first_derivative_fwd_y(
    Real top, Real middle, Real dy)
{
    return (top - middle) / dy;
}

inline Real first_derivative_bkwd_x(
    Real middle, Real left, Real dx)
{
    return (middle - left) / dx;
}

inline Real first_derivative_bkwd_y(
    Real middle, Real bottom, Real dy)
{
    return (middle - bottom) / dy;
}

inline Real first_derivative_of_square_x(
    Real right, Real middle, Real left, Real dx,
    Real alpha = 0.9)
{
    return 1./dx * (std::pow((middle + right) / 2., 2)
                    - std::pow((left + middle) / 2., 2))
            + alpha / dx
            * (std::abs(middle + right) * (middle - right) / 4.
                - std::abs(left + middle) * (left - middle) / 4.);
}

inline Real first_derivative_of_square_y(
    Real top, Real middle, Real bottom, Real dy,
    Real alpha = 0.9)
{
    return 1./dy * (std::pow((middle + top) / 2., 2)
                        - std::pow((bottom + middle) / 2., 2))
            + alpha / dy
            * (std::abs(middle + top) * (middle - top) / 4.
                - std::abs(bottom + middle) * (bottom - middle) / 4.);
}

inline Real first_derivative_of_product_x(
    Real u_left, Real u_middle, Real u_top, Real u_topleft,
    Real v_left, Real v_middle, Real v_right, Real dx,
    Real alpha = 0.9)
{
    return 1./dx * ((u_middle + u_top)
                        * (v_middle + v_right) / 4.
                        - (u_left + u_topleft) * (v_left + v_middle) / 4.)
            + alpha / dx * (std::abs(u_middle + u_top)
                * (v_middle - v_right) / 4.
            - std::abs(u_left + u_topleft) * (v_left - v_middle) / 4.);
}

inline Real first_derivative_of_product_y(
    Real v_right, Real v_middle, Real v_bottom,
    Real v_bottomright, Real u_bottom, Real u_middle,
    Real u_top, Real dy, Real alpha = 0.9)
{
    return 1./dy * ((v_middle + v_right)  * (u_middle + u_top) / 4. 
            - (v_bottom + v_bottomright) * (u_bottom + u_middle) / 4.)
            + alpha / dy * (std::abs(v_middle + v_right) 
                * (u_middle - u_top) / 4.
            - std::abs(v_bottom + v_bottomright) * (u_bottom - u_middle) / 4.);
}

inline Real first_derivative_u_temp_x(
    Real u_center, Real u_left, Real t_center, Real t_right,
    Real t_left, Real dx, Real alpha = 0.9)
{
    return 1./dx * (u_center * (t_center + t_right) / 2. 
                - u_left * (t_left + t_center) / 2.)
            + alpha / dx * (std::abs(u_center)*(t_center - t_right) / 2.
            - std::abs(u_left) * (t_left - t_center) / 2.);
}

inline Real first_derivative_v_temp_y(
    Real v_center, Real v_bottom, Real t_center, Real t_top, 
    Real t_bottom, Real dy, Real alpha = 0.9)
{
    return 1./dy * (v_center * (t_center + t_top)/ 2.
                - v_bottom * (t_bottom + t_center)/2.)
            + alpha / dy
                * (std::abs(v_center)*(t_center - t_top)/2.
                - std::abs(v_bottom)*(t_bottom - t_center)/2.);
}

}
}//namespace computation
#endif
