#ifndef NAST_HPX_GRID_STENCILS_HPP
#define NAST_HPX_GRID_STENCILS_HPP

#include <vector>

#include "partition_data.hpp"
#include "util/fd_stencils.hpp"
#include "util/cancellation_token.hpp"

namespace nast_hpx { namespace grid {

    static const std::size_t STENCIL_NONE = 20;
    static const std::size_t STENCIL_SET_VELOCITY_BOUNDARY = 21;
    static const std::size_t STENCIL_SET_VELOCITY_OBSTACLE = 22;
    static const std::size_t STENCIL_SET_VELOCITY_FLUID = 23;
    static const std::size_t STENCIL_COMPUTE_FG_BOUNDARY_AND_OBSTACLE = 24;
    static const std::size_t STENCIL_COMPUTE_FG_FLUID = 25;
    static const std::size_t STENCIL_COMPUTE_RHS = 26;
    static const std::size_t STENCIL_SET_P = 27;
    static const std::size_t STENCIL_SOR = 28;
    static const std::size_t STENCIL_JACOBI = 29;
    static const std::size_t STENCIL_COMPUTE_RESIDUAL = 30;
    static const std::size_t STENCIL_UPDATE_VELOCITY = 31;

    typedef std::pair<std::size_t, std::size_t> range_type;

    template <std::size_t stencil>
    struct stencils;

    template<>
    struct stencils<STENCIL_NONE>
    {
        static void call(partition_data<Real>& dst, partition_data<Real> const& src,
            std::vector<std::pair<std::size_t, std::size_t> > const& indices)
            {
                for (auto const& index_pair : indices)
                {
                    dst(index_pair.first, index_pair.second) =
                        src(index_pair.first, index_pair.second);
                }
            }
    };

    template<>
    struct stencils<STENCIL_SET_VELOCITY_BOUNDARY>
    {
        static void call(partition_data<Real>& dst_u, partition_data<Real>& dst_v,
            partition_data<Real> const& src_u, partition_data<Real> const& src_v,
            partition_data<std::bitset<6> > const& cell_types,
            std::vector<std::pair<std::size_t, std::size_t> > boundary_cells,
            boundary_data u_bnd, boundary_data v_bnd, boundary_type bnd_type)
            {
                for (auto const& idx_pair : boundary_cells)
                {
                    auto x = idx_pair.first;
                    auto y = idx_pair.second;

                    auto const& cell_type = cell_types(x, y);

                    if (cell_type.test(has_fluid_east))
                    {
                        switch(bnd_type.left)
                        {
                        case noslip:
                                dst_u(x, y) = 0;
                                dst_v(x, y) = -src_v(x + 1, y);
                                break;
                        case slip:
                                dst_u(x, y) = 0;
                                dst_v(x, y) = src_v(x + 1, y);
                                break;
                        case outstream:
                                dst_u(x, y) = src_u(x + 1, y);
                                dst_v(x, y) = src_v(x + 1, y);
                                break;
                        case instream:
                                dst_u(x, y) = u_bnd.left;
                                dst_v(x, y) = 2 * v_bnd.left - src_v(x, y);
                                break;
                        }
                    }

                    else if (cell_type.test(has_fluid_west))
                    {
                        switch(bnd_type.right)
                        {
                        case noslip:
                                dst_u(x - 1, y) = 0;
                                dst_v(x, y) = -src_v(x - 1, y);
                                break;
                        case slip:
                                dst_u(x - 1, y) = 0;
                                dst_v(x, y) = src_v(x - 1, y);
                                break;
                        case outstream:
                                dst_u(x - 1, y) = dst_u(x - 2, y);
                                dst_v(x, y) = src_v(x - 1, y);
                                break;
                        case instream:
                                dst_u(x - 1, y) = u_bnd.right;
                                dst_v(x, y) = 2 * v_bnd.right - src_v(x - 1, y);
                                break;
                        }
                    }

                    else if (cell_type.test(has_fluid_north))
                    {
                        switch(bnd_type.bottom)
                        {
                        case noslip:
                                dst_u(x, y) = -src_u(x, y + 1);
                                dst_v(x, y) = 0;
                                break;
                        case slip:
                                dst_u(x, y) = src_u(x, y + 1);
                                dst_v(x, y) = 0;
                                break;
                        case outstream:
                                dst_u(x, y) = src_u(x, y + 1);
                                dst_v(x, y) = src_v(x, y + 1);
                                break;
                        case instream:
                                dst_u(x, y) = 2 * u_bnd.bottom - src_u(x, y + 1);
                                dst_v(x, y) = v_bnd.bottom;
                                break;
                        }
                    }

                    else if (cell_type.test(has_fluid_south))
                    {
                        switch(bnd_type.top)
                        {
                        case noslip:
                                dst_u(x, y) = 2 * u_bnd.top - src_u(x, y - 1);
                                dst_v(x, y - 1) = 0;
                                break;
                        case slip:
                                dst_u(x, y) = src_u(x, y - 1);
                                dst_v(x, y - 1) = 0;
                                break;
                        case outstream:
                                dst_u(x, y) = src_u(x, y - 1);
                                dst_v(x, y - 1) = src_v(x, y - 1);
                                break;
                        case instream:
                                dst_u(x, y) = 2 * u_bnd.top - src_u(x, y - 1);
                                dst_v(x, y - 1) = v_bnd.top;
                                break;
                        }
                    }
                }
            }
    };

    template<>
    struct stencils<STENCIL_SET_VELOCITY_OBSTACLE>
    {
        static void call(partition_data<Real>& dst_u, partition_data<Real>& dst_v,
            partition_data<Real> const& src_u, partition_data<Real> const& src_v,
            partition_data<std::bitset<6> > const& cell_types,
            std::vector<std::pair<std::size_t, std::size_t> > const& obstacle_cells)
            {
                for (auto const& idx_pair : obstacle_cells)
                {
                    auto const x = idx_pair.first;
                    auto const y = idx_pair.second;

                    auto const& cell_type = cell_types(x, y);

                    dst_u(x, y) =
                        - src_u(x, y - 1) * cell_type.test(has_fluid_south)
                        - src_u(x, y + 1) * cell_type.test(has_fluid_north);

                    dst_v(x, y) =
                        - src_v(x - 1, y) * cell_type.test(has_fluid_west)
                        - src_v(x + 1, y) * cell_type.test(has_fluid_east);

                }
            }
    };

    template<>
    struct stencils<STENCIL_SET_VELOCITY_FLUID>
    {
        static void call(partition_data<Real>& dst_u, partition_data<Real>& dst_v,
            partition_data<Real> const& src_u, partition_data<Real> const& src_v,
            partition_data<std::bitset<6> > const& cell_types,
            std::vector<std::pair<std::size_t, std::size_t> > const& fluid_cells)
            {
                for (auto const& idx_pair : fluid_cells)
                {
                    auto const x = idx_pair.first;
                    auto const y = idx_pair.second;

                    auto const& cell_type = cell_types(x, y);

                    // only set fluid cells adjacent to obstacle
                    // fluid adjacent to boundary is set above
                    if (!cell_type.test(has_fluid_east) && !cell_types(x + 1, y).test(is_boundary))
                        dst_u(x, y) = 0;

                    if (!cell_type.test(has_fluid_north) && !cell_types(x, y + 1).test(is_boundary))
                        dst_v(x, y) = 0;
                }
            }
    };

    template<>
    struct stencils<STENCIL_COMPUTE_FG_BOUNDARY_AND_OBSTACLE>
    {
        static void call(partition_data<Real>& dst_f, partition_data<Real>& dst_g,
            partition_data<Real> const& src_u, partition_data<Real> const& src_v,
            partition_data<std::bitset<6> > const& cell_types,
            std::vector<std::pair<std::size_t, std::size_t> > const& boundary_cells,
            std::vector<std::pair<std::size_t, std::size_t> > const& obstacle_cells
           )
        {
            for (auto const& idx_pair : boundary_cells)
            {
                auto const x = idx_pair.first;
                auto const y = idx_pair.second;
                
                auto const& cell_type = cell_types(x, y);

                if (cell_type.test(has_fluid_east))
                    dst_f(x, y) = src_u(x, y);

                if (cell_type.test(has_fluid_north))
                    dst_g(x, y) = src_v(x, y);
            }            
            
            for (auto const& idx_pair : obstacle_cells)
            {
                auto const x = idx_pair.first;
                auto const y = idx_pair.second;
                
                auto const& cell_type = cell_types(x, y);

                if (cell_type.test(has_fluid_east))
                    dst_f(x, y) = src_u(x, y);

                if (cell_type.test(has_fluid_north))
                    dst_g(x, y) = src_v(x, y);
            }
        }
    };
    
    template<>
    struct stencils<STENCIL_COMPUTE_FG_FLUID>
    {
        static void call(partition_data<Real>& dst_f, partition_data<Real>& dst_g,
            partition_data<Real> const& src_u, partition_data<Real> const& src_v,
            partition_data<std::bitset<6> > const& cell_types,
            std::vector<std::pair<std::size_t, std::size_t> > const& fluid_cells,
            Real re, Real gx, Real gy, Real beta, Real dx, Real dy, Real dt,
            Real alpha)
        {
            for (auto const& idx_pair : fluid_cells)
            {
                auto const x = idx_pair.first;
                auto const y = idx_pair.second;
                
                auto const& cell_type = cell_types(x, y);

                    
                dst_f(x, y) = src_u(x, y) + cell_type.test(has_fluid_east)
                    * ( dt * (
                        1. / re * (
                            util::fd_stencils::second_derivative_fwd_bkwd_x(
                                src_u(x + 1, y), src_u(x, y), src_u(x - 1, y), dx
                            )
                            + util::fd_stencils::second_derivative_fwd_bkwd_y(
                                src_u(x, y + 1), src_u(x, y), src_u(x, y - 1), dy
                            )
                        )

                        - util::fd_stencils::first_derivative_of_square_x(
                            src_u(x + 1, y), src_u(x, y), src_u(x - 1, y),
                            dx, alpha
                        )
                        - util::fd_stencils::first_derivative_of_product_y(
                            src_v(x + 1, y), src_v(x, y), src_v(x, y - 1),
                            src_v(x + 1, y - 1), src_u(x, y - 1),
                            src_u(x, y), src_u(x, y + 1), dy, alpha
                        )

                        + gx
                    )
               /* - beta * dt / 2.
                                  * (middle_temperature
                                        + right_temperature )
                                  * gx*/
                );

                dst_g(x, y) = src_v(x, y) + cell_type.test(has_fluid_north)
                    * ( dt * (
                        1. / re * (
                            util::fd_stencils::second_derivative_fwd_bkwd_x(
                                src_v(x + 1, y), src_v(x, y), src_v(x - 1, y), dx
                            )
                            + util::fd_stencils::second_derivative_fwd_bkwd_y(
                                src_v(x, y + 1), src_v(x, y), src_v(x, y - 1), dy
                            )
                        )

                        - util::fd_stencils::first_derivative_of_product_x(
                            src_u(x - 1, y), src_u(x, y), src_u(x, y + 1),
                            src_u(x - 1, y + 1), src_v(x - 1, y),
                            src_v(x, y), src_v(x + 1, y), dx, alpha
                        )
                        - util::fd_stencils::first_derivative_of_square_y(
                            src_v(x, y + 1), src_v(x, y), src_v(x, y - 1),
                            dy, alpha
                        )

                        + gy
                    )
               /* - beta * dt / 2.
                                        * (middle_temperature
                                            + top_temperature )
                                        * gy*/
                );
            }
        }
    };

    template<>
    struct stencils<STENCIL_COMPUTE_RHS>
    {
        static void call(partition_data<Real>& dst_rhs,
            partition_data<Real> const& src_f, partition_data<Real> const& src_g,
            partition_data<std::bitset<6> > const& cell_types,
            std::vector<std::pair<std::size_t, std::size_t> > const& fluid_cells,
            Real dx, Real dy, Real dt)
        {
            for (auto const& idx_pair : fluid_cells)
            {
                auto const x = idx_pair.first;
                auto const y = idx_pair.second;

                dst_rhs(x, y) =
                    1. / dt *   (
                                    (src_f(x, y) - src_f(x - 1, y)) / dx
                                    +
                                    (src_g(x, y) - src_g(x, y - 1)) / dy
                                );
            }
        }
    };

    template<>
    struct stencils<STENCIL_SET_P>
    {
        static void call(partition_data<Real>& dst_p,
            partition_data<std::bitset<6> > const& cell_types,
            std::vector<std::pair<std::size_t, std::size_t> > const& boundary_cells,
            std::vector<std::pair<std::size_t, std::size_t> > const& obstacle_cells,
            util::cancellation_token token)
        {
            if (!token.was_cancelled())
            {
                for (auto& idx_pair : boundary_cells)
                {
                    auto x = idx_pair.first;
                    auto y = idx_pair.second;

                    auto const& cell_type = cell_types(x, y);

                    dst_p(x, y) =
                        (
                            dst_p(x - 1, y) * cell_type.test(has_fluid_west)
                            + dst_p(x + 1, y) * cell_type.test(has_fluid_east)
                            + dst_p(x, y - 1) * cell_type.test(has_fluid_south)
                            + dst_p(x, y + 1) * cell_type.test(has_fluid_north)
                        )
                        /
                        (
                            cell_type.test(has_fluid_west)
                            + cell_type.test(has_fluid_east)
                            + cell_type.test(has_fluid_south)
                            + cell_type.test(has_fluid_north)
                        );
                }

                for (auto& idx_pair : obstacle_cells)
                {
                    auto x = idx_pair.first;
                    auto y = idx_pair.second;

                    auto const& cell_type = cell_types(x, y);

                    dst_p(x, y) =
                        (
                            dst_p(x - 1, y) * cell_type.test(has_fluid_west)
                            + dst_p(x + 1, y) * cell_type.test(has_fluid_east)
                            + dst_p(x, y - 1) * cell_type.test(has_fluid_south)
                            + dst_p(x, y + 1) * cell_type.test(has_fluid_north)
                        )
                        /
                        (
                            cell_type.test(has_fluid_west)
                            + cell_type.test(has_fluid_east)
                            + cell_type.test(has_fluid_south)
                            + cell_type.test(has_fluid_north)
                        );
                }
            }
        }
    };

    template<>
    struct stencils<STENCIL_SOR>
    {
        //TODO remove idx, idy (was for debug)
        static void call(partition_data<Real>& dst_p,
            partition_data<Real> const& src_rhs,
            std::vector<std::pair<std::size_t, std::size_t> > const& fluid_cells,
            Real part1, Real part2, Real dx_sq, Real dy_sq,
            util::cancellation_token token)
        {
            if (!token.was_cancelled())
                for (auto& idx_pair : fluid_cells)
                {
                    auto x = idx_pair.first;
                    auto y = idx_pair.second;

                    dst_p(x, y) =
                        part1 * dst_p(x, y)
                        + part2 * (
                                (dst_p(x + 1, y) + dst_p(x - 1, y)) / dx_sq
                                + (dst_p(x, y + 1) + dst_p(x, y - 1)) / dy_sq
                                - src_rhs(x, y)
                        );
                }
        }
    };

    template<>
    struct stencils<STENCIL_JACOBI>
    {
        //TODO remove idx, idy (was for debug)
        static void call(partition_data<Real>& dst_p,
            partition_data<Real> const& src_rhs,
            std::vector<std::pair<std::size_t, std::size_t> > const& fluid_cells,
            Real dx_sq, Real dy_sq, util::cancellation_token token)
        {
            if (!token.was_cancelled())
                for (auto& idx_pair : fluid_cells)
                {
                    auto x = idx_pair.first;
                    auto y = idx_pair.second;

                    dst_p(x, y) =
                         ( (dst_p(x + 1, y) + dst_p(x - 1, y)) * dy_sq
                            + (dst_p(x, y + 1) + dst_p(x, y - 1)) * dx_sq
                            - dx_sq * dy_sq * src_rhs(x, y))
                        /
                        (2 * (dx_sq + dy_sq));
                }
        }
    };

    template<>
    struct stencils<STENCIL_COMPUTE_RESIDUAL>
    {
        static Real call(partition_data<Real> const& src_p,
            partition_data<Real> const& src_rhs,
            std::vector<std::pair<std::size_t, std::size_t> > const& fluid_cells,
            Real over_dx_sq, Real over_dy_sq, util::cancellation_token token)
        {
            Real local_residual = 0;

            if (!token.was_cancelled())
                for (auto& idx_pair : fluid_cells)
                {
                    auto x = idx_pair.first;
                    auto y = idx_pair.second;

                    Real tmp =
                        (src_p(x + 1, y) - 2 * src_p(x, y) + src_p(x - 1, y)) * over_dx_sq
                        + (src_p(x, y + 1) - 2 * src_p(x, y) + src_p(x, y - 1)) * over_dy_sq
                        - src_rhs(x, y);

                    local_residual += std::pow(tmp, 2);
                }

            return local_residual;
        }
    };

    template<>
    struct stencils<STENCIL_UPDATE_VELOCITY>
    {
        static std::pair<Real, Real>  call(partition_data<Real>& dst_u,
            partition_data<Real>& dst_v,
            partition_data<Real> const& src_f,
            partition_data<Real> const& src_g,
            partition_data<Real> const& src_p,
            partition_data<std::bitset<6> > const& cell_types,
            std::vector<std::pair<std::size_t, std::size_t> > const& fluid_cells,
            Real dt, Real over_dx, Real over_dy
            )
        {
            Real max_u = 0;
            Real max_v = 0;

            for (auto const& idx_pair : fluid_cells)
            {
                auto const x = idx_pair.first;
                auto const y = idx_pair.second;
                
                auto const& cell_type = cell_types(x, y);
                
                if (cell_type.test(has_fluid_east))
                {
                    dst_u(x, y) = src_f(x, y) - dt * over_dx *
                        (src_p(x + 1, y) - src_p(x, y));

                    max_u = std::abs(dst_u(x, y)) > max_u ? std::abs(dst_u(x, y)) : max_u;
                }

                if (cell_type.test(has_fluid_north))
                {
                    dst_v(x, y) = src_g(x, y) - dt * over_dy *
                        (src_p(x, y + 1) - src_p(x, y));

                    max_v = std::abs(dst_v(x, y)) > max_v ? std::abs(dst_v(x, y)) : max_v;
                }
            }                   

            return std::make_pair(max_u, max_v);
        }
    };
}
}

#endif
