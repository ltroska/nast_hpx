/** All stencils are provided in this header as inline methods
 */
#ifndef NAST_HPX_UTIL_FD_STENCILS_HPP
#define NAST_HPX_UTIL_FD_STENCILS_HPP

#include <cmath>
#include <cstdlib>

#include "defines.hpp"
#include "../grid/partition_data.hpp"

namespace nast_hpx { namespace util { namespace fd_stencils {

inline Real second_derivative_fwd_bkwd_x(grid::partition_data<Real> u_grid, std::size_t i, std::size_t j, Real over_dx_sq)
{
    return (u_grid(i + 1, j) - 2 * u_grid(i, j) + u_grid(i - 1, j)) * over_dx_sq;
}

inline Real second_derivative_fwd_bkwd_y(grid::partition_data<Real> v_grid, std::size_t i, std::size_t j, Real over_dy_sq)
{
    return (v_grid(i, j + 1) - 2 * v_grid(i, j) + v_grid(i, j - 1)) * over_dy_sq;
}

inline Real first_derivative_of_square_x(grid::partition_data<Real> u_grid, std::size_t i, std::size_t j, Real over_dx, Real alpha = 0.9)
{
    return (
                over_dx * (
                            std::pow(u_grid(i, j) + u_grid(i + 1, j), 2)
                            - std::pow(u_grid(i - 1, j) + u_grid(i, j), 2)
                           )

                + alpha * over_dx * (
                                            std::abs( u_grid(i, j) + u_grid(i + 1, j) )
                                            * ( u_grid(i, j) - u_grid(i + 1, j) )
                                        -
                                            std::abs( u_grid(i - 1, j) + u_grid(i, j))
                                            * ( u_grid(i - 1, j) - u_grid(i, j) )
                                    )
                ) / 4.;
}

inline Real first_derivative_of_square_y(grid::partition_data<Real> v_grid, std::size_t i, std::size_t j, Real over_dy, Real alpha = 0.9)
{
    return (
                over_dy *   (
                                std::pow(v_grid(i, j) + v_grid(i, j + 1), 2)
                                - std::pow(v_grid(i, j - 1) + v_grid(i, j), 2)
                            )

                + alpha * over_dy * (
                                            std::abs( v_grid(i, j) + v_grid(i, j + 1) )
                                            * ( v_grid(i, j) - v_grid(i, j + 1) )
                                        -
                                            std::abs( v_grid(i, j - 1) + v_grid(i, j))
                                            * ( v_grid(i, j - 1) - v_grid(i, j) )
                                    )
            ) / 4.;
}

inline Real first_derivative_of_product_x(grid::partition_data<Real> u_grid, grid::partition_data<Real> v_grid, std::size_t i, std::size_t j, Real over_dx, Real alpha = 0.9)
{
    return over_dx * (      ( u_grid(i, j) + u_grid(i, j + 1))
                            * (v_grid(i, j) + v_grid(i + 1, j))
                            -
                                (u_grid(i - 1, j) + u_grid(i - 1, j + 1))
                                 * (v_grid(i - 1, j) + v_grid(i, j))
                      ) / 4.

           + alpha * over_dx * (
                                    std::abs(u_grid(i, j) + u_grid(i, j + 1))
                                    * (v_grid(i, j) - v_grid(i + 1, j))
                                    -
                                        std::abs(u_grid(i - 1, j) + u_grid(i - 1, j + 1))
                                        * (v_grid(i - 1, j) - v_grid(i, j))
                                ) / 4.;
}

inline Real first_derivative_of_product_y(grid::partition_data<Real> u_grid, grid::partition_data<Real> v_grid, std::size_t i, std::size_t j, Real over_dy, Real alpha = 0.9)
{
    return over_dy * (      ( v_grid(i, j) + v_grid(i + 1, j))
                            * (u_grid(i, j) + u_grid(i, j + 1))
                            -
                                (v_grid(i, j - 1) + v_grid(i + 1, j - 1))
                                 * (u_grid(i, j - 1) + u_grid(i, j))
                      ) / 4.

           + alpha * over_dy * (
                                    std::abs(v_grid(i, j) + v_grid(i + 1, j))
                                    * (u_grid(i, j) - u_grid(i, j + 1))
                                    -
                                        std::abs(v_grid(i, j - 1) + u_grid(i + 1, j - 1))
                                        * (u_grid(i, j - 1) - u_grid(i, j))
                                ) / 4.;
}

inline Real interpolate(Real x, Real y, Real x1, Real x2, Real y1, Real y2, Real v1, Real v2, Real v3, Real v4, Real dx, Real dy)
{
    return 1./(dx * dy) * ( (x2 - x)*(y2 - y) * v1 + (x - x1)*(y2 - y) * v2
                            + (x2 - x)*(y - y1) * v3 + (x - x1)*(y - y1) * v4);
}

}
}
}
#endif
