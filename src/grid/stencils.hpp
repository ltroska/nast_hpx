#ifndef NAST_HPX_GRID_STENCILS_HPP
#define NAST_HPX_GRID_STENCILS_HPP

#include <vector>

#include "partition_data.hpp"
#include "../util/fd_stencils.hpp"
#include "../util/cancellation_token.hpp"
#include "boundary_condition.hpp"

namespace nast_hpx { namespace grid {

    static const std::size_t STENCIL_NONE = 20;
    static const std::size_t STENCIL_SET_VELOCITY = 23;
    static const std::size_t STENCIL_UPDATE_PARTICLE_POSITION = 24;
    static const std::size_t STENCIL_COMPUTE_FG = 33;
    static const std::size_t STENCIL_COMPUTE_RHS = 26;
    static const std::size_t STENCIL_SET_P = 27;
    static const std::size_t STENCIL_SOR = 28;
    static const std::size_t STENCIL_JACOBI = 29;
    static const std::size_t STENCIL_COMPUTE_RESIDUAL = 30;
    static const std::size_t STENCIL_UPDATE_VELOCITY = 31;

    typedef index range_type;

    template <std::size_t stencil>
    struct stencils;

    template<>
    struct stencils<STENCIL_NONE>
    {
        static void call(partition_data<Real>& dst, partition_data<Real> const& src,
            std::vector<index> const& indices)
            {
                for (auto const& index_pair : indices)
                {
                    dst(index_pair.x, index_pair.y) =
                        src(index_pair.x, index_pair.y);
                }
            }
    };

    template<>
    struct stencils<STENCIL_SET_VELOCITY>
    {
        static void call(partition_data<Real>& dst_u, partition_data<Real>& dst_v,
            partition_data<std::bitset<7> > const& cell_types,
            std::vector<index> const& obstacle_cells,
            boundary_condition const& bnd_condition)
            {

                for (auto const& idx_pair : obstacle_cells)
                {
                    auto i = idx_pair.x;
                    auto j = idx_pair.y;

                    auto const& cell_type = cell_types(i, j);

                    if (cell_type.test(is_boundary))
                    {
                        if (cell_type.test(has_fluid_right))
                        {
                            switch(bnd_condition.left_type)
                            {
                            case noslip:
                                    dst_u(i, j) = 0;
                                    dst_v(i, j) = -dst_v(i + 1, j);
                                    break;
                            case slip:
                                    dst_u(i, j) = 0;
                                    dst_v(i, j) = dst_v(i + 1, j);
                                    break;
                            case outstream:
                                    dst_u(i, j) = dst_u(i + 1, j);
                                    dst_v(i, j) = dst_v(i + 1, j);
                                    break;
                            case instream:
                                    dst_u(i, j) = bnd_condition.left.x;
                                    dst_v(i, j) = 2 * bnd_condition.left.y - dst_v(i + 1, j);
                                    break;
                            }
                        }

                        else if (cell_type.test(has_fluid_left))
                        {
                            switch(bnd_condition.right_type)
                            {
                            case noslip:
                                    dst_u(i - 1, j) = 0;
                                    dst_v(i, j) = -dst_v(i - 1, j);
                                    break;
                            case slip:
                                    dst_u(i - 1, j) = 0;
                                    dst_v(i, j) = dst_v(i - 1, j);
                                    break;
                            case outstream:
                                    dst_u(i - 1, j) = dst_u(i - 2, j);
                                    dst_v(i, j) = dst_v(i - 1, j);
                                    break;
                            case instream:
                                    dst_u(i - 1, j) = bnd_condition.right.x;
                                    dst_v(i, j) = 2 * bnd_condition.right.y - dst_v(i - 1, j);
                                    break;
                            }
                        }

                        else if (cell_type.test(has_fluid_top))
                        {
                            switch(bnd_condition.bottom_type)
                            {
                            case noslip:
                                    dst_u(i, j) = -dst_u(i, j + 1);
                                    dst_v(i, j) = 0;
                                    break;
                            case slip:
                                    dst_u(i, j) = dst_u(i, j + 1);
                                    dst_v(i, j) = 0;
                                    break;
                            case outstream:
                                    dst_u(i, j) = dst_u(i, j + 1);
                                    dst_v(i, j) = dst_v(i, j + 1);
                                    break;
                            case instream:
                                    dst_u(i, j) = 2 * bnd_condition.bottom.x - dst_u(i, j + 1);
                                    dst_v(i, j) = bnd_condition.bottom.y;
                                    break;
                            }
                        }

                        else if (cell_type.test(has_fluid_bottom))
                        {
                            switch(bnd_condition.top_type)
                            {
                            case noslip:
                                    dst_u(i, j) = 2 * bnd_condition.top.x - dst_u(i, j - 1);
                                    dst_v(i, j - 1) = 0;
                                    break;
                            case slip:
                                    dst_u(i, j) = dst_u(i, j - 1);
                                    dst_v(i, j - 1) = 0;
                                    break;
                            case outstream:
                                    dst_u(i, j) = dst_u(i, j - 1);
                                    dst_v(i, j - 1) = dst_v(i, j - 2);
                                    break;
                            case instream:
                                    dst_u(i, j) = 2 * bnd_condition.top.x - dst_u(i, j - 1);
                                    dst_v(i, j - 1) = bnd_condition.top.y;
                                    break;
                            }
                        }
                    }
                    else
                    {
                        dst_u(i, j) =
                            - dst_u(i, j - 1) * cell_type.test(has_fluid_bottom)
                            - dst_u(i, j + 1) * cell_type.test(has_fluid_top);

                        dst_v(i, j) =
                            - dst_v(i - 1, j) * cell_type.test(has_fluid_left)
                            - dst_v(i + 1, j) * cell_type.test(has_fluid_right);
                    }
                }
            }
    };

    template<>
    struct stencils<STENCIL_UPDATE_PARTICLE_POSITION>
    {
        static void call(std::vector<pair<Real> >::iterator particlesBegin,
            std::vector<pair<Real> >::iterator particlesEnd,
            partition_data<Real> const& src_u, partition_data<Real> const& src_v,
            Real dx, Real dy, Real dt, Real x_length, Real y_length)
        {
            for (auto it = particlesBegin; it < particlesEnd; ++it)
            {
                auto const x_pos = it->x;
                auto const y_pos = it->y;

                std::size_t i = (int) (x_pos/dx) + 1;
                std::size_t j = (int) ((y_pos + dy/2.)/dy) + 1;

                Real x1 = (i - 1) * dx;
                Real x2 = i * dx;
                Real y1 = (j - 1.5) * dy;
                Real y2 = (j - 0.5) * dy;

                it->x = std::max(0., std::min(x_length, x_pos + dt * util::fd_stencils::interpolate(x_pos, y_pos, x1, x2, y1, y2,
                                                       src_u(i - 1, j - 1), src_u(i, j - 1), src_u(i - 1, j), src_u(i, j)
                                                       , dx, dy)));

                i = (int) ((x_pos + dx/2.)/dx) + 1;
                j = (int) (y_pos/dy) + 1;

                x1 = (i - 1.5) * dx;
                x2 = (i - 0.5) * dx;
                y1 = (j - 1) * dy;
                y2 = j * dy;

                it->y = std::max(0., std::min(y_length, y_pos + dt * util::fd_stencils::interpolate(x_pos, y_pos, x1, x2, y1, y2,
                                                      src_v(i - 1, j - 1), src_v(i, j - 1), src_v(i - 1, j), src_v(i, j)
                                                      , dx, dy)));
            }
        }
    };

    template<>
    struct stencils<STENCIL_COMPUTE_FG>
    {
        static void call(partition_data<Real>& dst_f, partition_data<Real>& dst_g,
            partition_data<Real> const& src_u, partition_data<Real> const& src_v,
            partition_data<std::bitset<7> > const& cell_types,
            std::vector<index> const& obstacle_cells,
            std::vector<index> const& fluid_cells,
            Real re, Real gx, Real gy, Real dx, Real dy, Real dt,
            Real alpha
           )
        {


            Real dx_sq = std::pow(dx, 2);
            Real dy_sq = std::pow(dx, 2);

            Real over_dx = 1./dx;
            Real over_dy = 1./dy;

            Real over_dx_sq = 1./dx_sq;
            Real over_dy_sq = 1./dy_sq;

            for (auto const& idx_pair : obstacle_cells)
            {
                auto const i = idx_pair.x;
                auto const j = idx_pair.y;

                auto const& cell_type = cell_types(i, j);

                if (cell_type.test(has_fluid_right))
                    dst_f(i, j) = src_u(i, j);

                if (cell_type.test(has_fluid_top))
                    dst_g(i, j) = src_v(i, j);
            }

            for (auto const& idx_pair : fluid_cells)
            {
                auto const i = idx_pair.x;
                auto const j = idx_pair.y;

                auto const& cell_type = cell_types(i, j);


                dst_f(i, j) = src_u(i, j) + cell_type.test(has_fluid_right)
                    * ( dt * (
                        1. / re * (
                            util::fd_stencils::second_derivative_fwd_bkwd_x(src_u, i, j, over_dx_sq)
                            + util::fd_stencils::second_derivative_fwd_bkwd_y(src_u, i, j, over_dy_sq)
                        )

                        - util::fd_stencils::first_derivative_of_square_x(src_u, i, j, over_dx, alpha)
                        - util::fd_stencils::first_derivative_of_product_y(src_u, src_v, i, j, over_dy, alpha)
                        + gx
                        )
                    );

                dst_g(i, j) = src_v(i, j) + cell_type.test(has_fluid_top)
                    * ( dt * (
                        1. / re * (
                            util::fd_stencils::second_derivative_fwd_bkwd_x(src_v, i, j, over_dx_sq)
                            + util::fd_stencils::second_derivative_fwd_bkwd_y(src_v, i, j, over_dy_sq)
                        )

                        - util::fd_stencils::first_derivative_of_product_x(src_u, src_v, i, j, over_dx, alpha)
                        - util::fd_stencils::first_derivative_of_square_y(src_v, i, j, over_dy, alpha)
                        + gy
                        )
                    );
            }
        }
    };


    template<>
    struct stencils<STENCIL_COMPUTE_RHS>
    {
        static void call(partition_data<Real>& dst_rhs,
            partition_data<Real> const& src_f, partition_data<Real> const& src_g,
            partition_data<std::bitset<7> > const& cell_types,
            std::vector<index> const& fluid_cells,
            Real dx, Real dy, Real dt)
        {


            for (auto const& idx_pair : fluid_cells)
            {
                auto const i = idx_pair.x;
                auto const j = idx_pair.y;

                dst_rhs(i, j) =
                    1. / dt *   (
                                    (src_f(i, j) - src_f(i - 1, j)) / dx
                                    +
                                    (src_g(i, j) - src_g(i, j - 1)) / dy
                                );
            }
        }
    };

    template<>
    struct stencils<STENCIL_SET_P>
    {
        static void call(partition_data<Real>& dst_p,
            partition_data<std::bitset<7> > const& cell_types,
            std::vector<index> const& obstacle_cells,
            util::cancellation_token token)
        {


            if (!token.was_cancelled())
            {
                for (auto& idx_pair : obstacle_cells)
                {
                    auto i = idx_pair.x;
                    auto j = idx_pair.y;

                    auto const& cell_type = cell_types(i, j);

                    dst_p(i, j) =
                        (
                            dst_p(i - 1, j) * cell_type.test(has_fluid_left)
                            + dst_p(i + 1, j) * cell_type.test(has_fluid_right)
                            + dst_p(i, j - 1) * cell_type.test(has_fluid_bottom)
                            + dst_p(i, j + 1) * cell_type.test(has_fluid_top)
                        )
                        /
                        (
                            cell_type.test(has_fluid_left)
                            + cell_type.test(has_fluid_right)
                            + cell_type.test(has_fluid_bottom)
                            + cell_type.test(has_fluid_top)
                        );
                }
            }
        }
    };

    template<>
    struct stencils<STENCIL_SOR>
    {
        static void call(partition_data<Real>& dst_p,
            partition_data<Real> const& src_rhs,
            std::vector<index> const& fluid_cells,
            Real part1, Real part2, Real dx_sq, Real dy_sq,
            util::cancellation_token token)
        {

            if (!token.was_cancelled())
            {
                for (auto& idx_pair : fluid_cells)
                {
                    auto i = idx_pair.x;
                    auto j = idx_pair.y;

                    dst_p(i, j) =
                        part1 * dst_p(i, j)
                        + part2 * (
                                (dst_p(i + 1, j) + dst_p(i - 1, j)) / dx_sq
                                + (dst_p(i, j + 1) + dst_p(i, j - 1)) / dy_sq
                                - src_rhs(i, j)
                        );
                }
            }
        }
    };

    template<>
    struct stencils<STENCIL_JACOBI>
    {
        static void call(partition_data<Real>& dst_p,
            partition_data<Real> const& src_rhs,
            std::vector<index> const& fluid_cells,
            Real dx_sq, Real dy_sq, util::cancellation_token token)
        {


            if (!token.was_cancelled())
            {
                for (auto& idx_pair : fluid_cells)
                {
                    auto i = idx_pair.x;
                    auto j = idx_pair.y;

                    dst_p(i, j) =
                         ( (dst_p(i + 1, j) + dst_p(i - 1, j)) * dy_sq
                            + (dst_p(i, j + 1) + dst_p(i, j - 1)) * dx_sq
                            - dx_sq * dy_sq * src_rhs(i, j))
                        /
                        (2 * (dx_sq + dy_sq));
                }
            }
        }
    };

    template<>
    struct stencils<STENCIL_COMPUTE_RESIDUAL>
    {
        static Real call(partition_data<Real> const& src_p,
            partition_data<Real> const& src_rhs,
            std::vector<index> const& fluid_cells,
            Real over_dx_sq, Real over_dy_sq, util::cancellation_token token)
        {

            Real local_residual = 0;

            if (!token.was_cancelled())
                for (auto& idx_pair : fluid_cells)
                    {
                        auto i = idx_pair.x;
                        auto j = idx_pair.y;

                        Real tmp =
                            (src_p(i + 1, j) - 2 * src_p(i, j) + src_p(i - 1, j)) * over_dx_sq
                            + (src_p(i, j + 1) - 2 * src_p(i, j) + src_p(i, j - 1)) * over_dy_sq
                            - src_rhs(i, j);

                        local_residual += std::pow(tmp, 2);
                    }

            return local_residual;
        }
    };

    template<>
    struct stencils<STENCIL_UPDATE_VELOCITY>
    {
        static pair<Real> call(partition_data<Real>& dst_u,
            partition_data<Real>& dst_v,
            partition_data<Real> const& src_f,
            partition_data<Real> const& src_g,
            partition_data<Real> const& src_p,
            partition_data<std::bitset<7> > const& cell_types,
            std::vector<index> const& fluid_cells,
            Real dt, Real over_dx, Real over_dy
            )
        {


            pair<Real> max_uv(0);

            for (auto const& idx_pair : fluid_cells)
            {
                auto const i = idx_pair.x;
                auto const j = idx_pair.y;

                auto const& cell_type = cell_types(i, j);

                if (cell_type.test(has_fluid_right))
                {
                    dst_u(i, j) = src_f(i, j) - dt * over_dx *
                        (src_p(i + 1, j) - src_p(i, j));

                    max_uv.x = std::abs(dst_u(i, j)) > max_uv.x ? std::abs(dst_u(i, j)) : max_uv.x;
                }

                if (cell_type.test(has_fluid_top))
                {
                    dst_v(i, j) = src_g(i, j) - dt * over_dy *
                        (src_p(i, j + 1) - src_p(i, j));

                    max_uv.y = std::abs(dst_v(i, j)) > max_uv.y ? std::abs(dst_v(i, j)) : max_uv.y;
                }
            }

            return max_uv;
        }
    };
}
}

#endif
