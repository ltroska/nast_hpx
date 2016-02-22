#ifndef COMPUTATION_STENCILS_HPP
#define COMPUTATION_STENCILS_HPP

#include "util/types.hpp"

namespace computation {

RealType second_derivative_fwd_bkwd_x(RealType right, RealType middle, RealType left, RealType dx)
{
    return (right-2*middle+left)/(dx*dx);
}

RealType second_derivative_fwd_bkwd_y(RealType top, RealType middle, RealType bottom, RealType dy)
{
    return (top-2*middle+bottom)/(dy*dy);
}

RealType first_derivative_fwd_x(RealType right, RealType middle, RealType dx)
{
    return (right-middle)/dx;
}

RealType first_derivative_fwd_y(RealType top, RealType middle, RealType dy)
{
    return (top-middle)/dy;
}

RealType first_derivative_bkwd_x(RealType middle, RealType left, RealType dx)
{
    return (middle-left)/dx;
}

RealType first_derivative_bkwd_y(RealType middle, RealType bottom, RealType dy)
{
    return (middle-bottom)/dy;
}

RealType first_derivative_of_square_x(RealType right, RealType middle, RealType left, RealType dx, RealType alpha = 0.9)
{
    return 1./dx*(((middle+right)/2.)*((middle+right)/2.)-((left+middle)/2.)*((left+middle)/2.))
            +alpha/dx*(abs(middle+right)/2 * (middle-right)/2 - abs(left+middle)/2 * (left-middle)/2);
}

RealType first_derivative_of_square_y(RealType top, RealType middle, RealType bottom, RealType dy, RealType alpha = 0.9)
{
    return 1./dy*(((middle+top)/2.)*((middle+top)/2.)-((bottom+middle)/2.)*((bottom+middle)/2.))
            +alpha/dy*(abs(middle+top)/2 * (middle-top)/2 - abs(bottom+middle)/2 * (bottom-middle)/2);
}

RealType first_derivative_of_product_x(RealType u_left, RealType u_middle, RealType u_top, RealType u_topleft,
                                       RealType v_left, RealType v_middle, RealType v_right, RealType dx, RealType alpha = 0.9)
{
    return 1./dx*( ((u_middle+u_top)/2)*((u_middle+u_left)/2) - ((u_left+u_topleft)/2)*((v_left+v_middle)/2) )
            +alpha/dx*(abs(u_middle+u_top)/2*(v_middle-v_right)/2 - abs(u_left+u_topleft)/2 * (v_left-v_middle)/2);
}

RealType first_derivative_of_product_y(RealType v_right, RealType v_middle, RealType v_bottom, RealType v_bottomright,
                                       RealType u_bottom, RealType u_middle, RealType u_top, RealType dy, RealType alpha = 0.9)
{
    return 1./dy*((v_middle+v_right)/2*(u_middle+u_top)/2 - (v_bottom+v_bottomright)/2*(u_bottom+u_middle)/2)
            +alpha/dy*(abs(v_middle+v_right)/2*(u_middle-u_top)/2 - abs(v_bottom+v_bottomright)/2*(u_bottom-u_middle));
}


}//namespace computation
#endif
