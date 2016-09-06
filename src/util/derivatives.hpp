/** All stencils are provided in this header as inline methods
 */
#ifndef NAST_HPX_UTIL_DERIVATIVES_HPP
#define NAST_HPX_UTIL_DERIVATIVES_HPP

#include <cmath>
#include <cstdlib>

#include "grid/partition_data.hpp"

namespace nast_hpx { namespace util { namespace derivatives {

typedef grid::partition_data<Real> grid_type;

inline Real second_derivative_fwd_bkwd_x(
    grid_type grid, std::size_t i, std::size_t j, std::size_t k, Real dx_sq)
{
    return (grid(i + 1, j, k) - 2 * grid(i, j, k) + grid(i - 1, j , k)) / dx_sq;
}

inline Real second_derivative_fwd_bkwd_y(
    grid_type grid, std::size_t i, std::size_t j, std::size_t k, Real dy_sq)
{
    return (grid(i, j + 1, k) - 2 * grid(i, j, k) + grid(i, j - 1, k)) / dy_sq;
}

inline Real second_derivative_fwd_bkwd_z(
    grid_type grid, std::size_t i, std::size_t j, std::size_t k, Real dz_sq)
{
    return (grid(i, j, k + 1) - 2 * grid(i, j, k) + grid(i, j, k - 1)) / dz_sq;
}

inline Real first_derivative_fwd_x(
    grid_type grid, std::size_t i, std::size_t j, std::size_t k, Real dx)
{
    return (grid(i + 1, j, k) - grid(i, j, k)) / dx;
}

inline Real first_derivative_fwd_y(
    grid_type grid, std::size_t i, std::size_t j, std::size_t k, Real dy)
{
    return (grid(i, j + 1, k) - grid(i, j, k)) / dy;
}

inline Real first_derivative_fwd_z(
    grid_type grid, std::size_t i, std::size_t j, std::size_t k, Real dz)
{
    return (grid(i, j, k + 1) - grid(i, j, k)) / dz;
}

inline Real first_derivative_bkwd_x(
    grid_type grid, std::size_t i, std::size_t j, std::size_t k, Real dx)
{
    return (grid(i, j, k) - grid(i - 1, j, k)) / dx;
}

inline Real first_derivative_bkwd_y(
    grid_type grid, std::size_t i, std::size_t j, std::size_t k, Real dy)
{
    return (grid(i, j, k) - grid(i, j - 1, k)) / dy;
}

inline Real first_derivative_bkwd_z(
    grid_type grid, std::size_t i, std::size_t j, std::size_t k, Real dz)
{
    return (grid(i, j, k) - grid(i, j, k - 1)) / dz;
}

inline Real first_derivative_of_square_x(
    grid_type grid, std::size_t i, std::size_t j, std::size_t k, Real dx,
    Real alpha = 0.9)
{
    return 1./dx * (std::pow((grid(i, j, k) + grid(i + 1, j, k)) / 2., 2)
                    - std::pow((grid(i - 1, j, k) + grid(i, j, k)) / 2., 2))
            + alpha / dx
            * (std::abs(grid(i, j, k) + grid(i + 1, j, k)) * (grid(i, j, k) - grid(i + 1, j, k)) / 4.
                - std::abs(grid(i - 1, j, k) + grid(i, j, k)) * (grid(i - 1, j, k) - grid(i, j, k)) / 4.);
}

inline Real first_derivative_of_square_y(
    grid_type grid, std::size_t i, std::size_t j, std::size_t k, Real dy,
    Real alpha = 0.9)
{
    return 1./dy * (std::pow((grid(i, j, k) + grid(i, j + 1, k)) / 2., 2)
                        - std::pow((grid(i, j - 1, k) + grid(i, j, k)) / 2., 2))
            + alpha / dy
            * (std::abs(grid(i, j, k) + grid(i, j + 1, k)) * (grid(i, j, k) - grid(i, j + 1, k)) / 4.
                - std::abs(grid(i, j - 1, k) + grid(i, j, k)) * (grid(i, j - 1, k) - grid(i, j, k)) / 4.);
}

inline Real first_derivative_of_square_z(
    grid_type grid, std::size_t i, std::size_t j, std::size_t k, Real dz,
    Real alpha = 0.9)
{
    return 1./dz * (std::pow((grid(i, j, k) + grid(i, j, k + 1)) / 2., 2)
                        - std::pow((grid(i, j, k - 1) + grid(i, j, k)) / 2., 2))
            + alpha / dz
            * (std::abs(grid(i, j, k) + grid(i, j, k + 1)) * (grid(i, j, k) - grid(i, j, k + 1)) / 4.
                - std::abs(grid(i, j, k - 1) + grid(i, j, k)) * (grid(i, j, k - 1) - grid(i, j, k)) / 4.);
}

inline Real first_derivative_of_uv_x(
    grid_type u, grid_type v, std::size_t i, std::size_t j, std::size_t k, Real dx,
    Real alpha = 0.9)
{
    return 1./dx * ((u(i, j, k) + u(i, j + 1, k))
                        * (v(i, j, k) + v(i + 1, j, k)) / 4.
                        - (u(i - 1, j, k) + u(i - 1, j + 1, k)) * (v(i - 1, j , k) + v(i, j, k)) / 4.)
            + alpha / dx * (std::abs(u(i, j, k) + u(i, j + 1, k))
                * (v(i, j, k) - v(i + 1, j, k)) / 4.
            - std::abs(u(i - 1, j, k) + u(i - 1, j + 1, k)) * (v(i - 1, j, k) - v(i, j, k)) / 4.);
}

inline Real first_derivative_of_uw_x(
    grid_type u, grid_type w, std::size_t i, std::size_t j, std::size_t k, Real dx,
    Real alpha = 0.9)
{
    return 1./dx * ((u(i, j, k) + u(i, j, k + 1))
                        * (w(i, j, k) + w(i + 1, j, k)) / 4.
                        - (u(i - 1, j, k) + u(i - 1, j, k + 1)) * (w(i - 1, j , k) + w(i, j, k)) / 4.)
            + alpha / dx * (std::abs(u(i, j, k) + u(i, j, k + 1))
                * (w(i, j, k) - w(i + 1, j, k)) / 4.
            - std::abs(u(i - 1, j, k) + u(i - 1, j, k + 1)) * (w(i - 1, j, k) - w(i, j, k)) / 4.);
}

inline Real first_derivative_of_uv_y(
    grid_type u, grid_type v, std::size_t i, std::size_t j, std::size_t k, Real dy
    , Real alpha = 0.9)
{
    return 1./dy * ((v(i, j, k) + v(i + 1, j, k))  * (u(i, j, k) + u(i, j + 1, k)) / 4.
            - (v(i, j - 1, k) + v(i + 1, j - 1, k)) * (u(i, j - 1, k) + u(i, j, k)) / 4.)
            + alpha / dy * (std::abs(v(i, j, k) + v(i + 1, j, k))
                * (u(i, j, k) - u(i, j + 1, k)) / 4.
            - std::abs(v(i, j - 1, k) + v(i + 1, j - 1, k)) * (u(i, j - 1, k) - u(i, j, k)) / 4.);
}

inline Real first_derivative_of_vw_y(
    grid_type v, grid_type w, std::size_t i, std::size_t j, std::size_t k, Real dy
    , Real alpha = 0.9)
{
    return 1./dy * ((v(i, j, k) + v(i, j, k + 1))  * (w(i, j, k) + w(i, j + 1, k)) / 4.
            - (v(i, j - 1, k) + v(i, j - 1, k + 1)) * (w(i, j - 1, k) + w(i, j, k)) / 4.)
            + alpha / dy * (std::abs(v(i, j, k) + v(i, j, k + 1))
                * (w(i, j, k) - w(i, j + 1, k)) / 4.
            - std::abs(v(i, j - 1, k) + v(i, j - 1, k + 1)) * (w(i, j - 1, k) - w(i, j, k)) / 4.);
}

inline Real first_derivative_of_uw_z(
    grid_type u, grid_type w, std::size_t i, std::size_t j, std::size_t k, Real dz
    , Real alpha = 0.9)
{
    return 1./dz * ((w(i, j, k) + w(i + 1, j, k))  * (u(i, j, k) + u(i, j, k + 1)) / 4.
            -
            (w(i, j, k - 1) + w(i + 1, j, k - 1)) * (u(i, j, k - 1) + u(i, j, k)) / 4.)
            +
            alpha / dz *
                (std::abs(w(i, j, k) + w(i + 1, j, k)) * (u(i, j, k) - u(i, j, k + 1)) / 4.
            - std::abs(w(i, j, k - 1) + w(i + 1, j, k - 1)) * (u(i, j, k - 1) - u(i, j, k)) / 4.);
}

inline Real first_derivative_of_vw_z(
    grid_type v, grid_type w, std::size_t i, std::size_t j, std::size_t k, Real dz
    , Real alpha = 0.9)
{
    return 1./dz * ((w(i, j, k) + w(i, j + 1, k))  * (v(i, j, k) + v(i, j, k + 1)) / 4.
            - (w(i, j, k - 1) + w(i, j + 1, k - 1)) * (v(i, j, k - 1) + v(i, j, k)) / 4.)
            + alpha / dz * (std::abs(w(i, j, k) + w(i, j + 1, k))
                * (v(i, j, k) - v(i, j, k + 1)) / 4.
            - std::abs(w(i, j, k - 1) + w(i, j + 1, k - 1)) * (v(i, j, k - 1) - v(i, j, k)) / 4.);
}

}
}
}
#endif
