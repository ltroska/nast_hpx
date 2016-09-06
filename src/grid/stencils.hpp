#ifndef NAST_HPX_GRID_STENCILS_HPP
#define NAST_HPX_GRID_STENCILS_HPP

#include <vector>

#ifdef WITH_FOR_EACH
#include <hpx/parallel/algorithms/for_each.hpp>
#include <hpx/parallel/algorithms/reduce.hpp>
#include <hpx/parallel/algorithms/transform_reduce.hpp>
#endif

#include "partition_data.hpp"
#include "util/derivatives.hpp"
#include "util/cancellation_token.hpp"

namespace nast_hpx { namespace grid {

    static const std::size_t STENCIL_NONE = 20;
    static const std::size_t STENCIL_SET_VELOCITY_OBSTACLE = 23;
    static const std::size_t STENCIL_SET_FG_OBSTACLE = 24;
    static const std::size_t STENCIL_COMPUTE_FG = 33;
    static const std::size_t STENCIL_COMPUTE_RHS = 26;
    static const std::size_t STENCIL_SET_P_OBSTACLE = 27;
    static const std::size_t STENCIL_SOR = 28;
    static const std::size_t STENCIL_JACOBI = 29;
    static const std::size_t STENCIL_COMPUTE_RESIDUAL = 30;
    static const std::size_t STENCIL_UPDATE_VELOCITY = 31;
    static const std::size_t STENCIL_TEST = 32;

    typedef std::pair<std::size_t, std::size_t> range_type;

    template <std::size_t stencil>
    struct stencils;

    template<>
    struct stencils<STENCIL_NONE>
    {
        static void call(partition_data<Real>& dst, partition_data<Real> const& src,
                std::vector<index> const& indices)
            {
                for (auto const& ind : indices)
                {
                    dst(ind.x, ind.y, ind.z) =
                        src(ind.x, ind.y, ind.z);
                }
            }
    };


    template<>
    struct stencils<STENCIL_SET_VELOCITY_OBSTACLE>
    {
        static void call(partition_data<Real>& dst_u, partition_data<Real>& dst_v,
                partition_data<Real>& dst_w,
                partition_data<std::bitset<9> > const& cell_types,
                std::vector<index> const& boundary_cells,
                std::vector<index> const& obstacle_cells,
                boundary_condition const& bnd_condition)
            {
                for  (auto const& ind : boundary_cells)
                {
                    auto const i = ind.x;
                    auto const j = ind.y;
                    auto const k = ind.z;

                    auto const& cell_type = cell_types(i, j, k);

                    if (cell_type.test(has_fluid_bottom))
                    {
                        switch(bnd_condition.top_type)
                        {
                            case noslip:
                                dst_u(i, j, k) = 2 * bnd_condition.top.x - dst_u(i, j, k - 1);
                                dst_v(i, j, k) = 2 * bnd_condition.top.y - dst_v(i, j, k - 1);
                                dst_w(i, j, k - 1) = 0;
                                break;
                        }
                    }

                    if (cell_type.test(has_fluid_top))
                    {
                        switch(bnd_condition.bottom_type)
                        {
                            case noslip:
                                dst_u(i, j, k) = 2 * bnd_condition.bottom.x - dst_u(i, j, k + 1);
                                dst_v(i, j, k) = 2 * bnd_condition.bottom.y - dst_v(i, j, k + 1);
                                dst_w(i, j, k) = 0;
                                break;

                        }
                    }

                    if (cell_type.test(has_fluid_left))
                    {
                        switch(bnd_condition.right_type)
                        {
                            case noslip:
                                dst_u(i - 1, j, k) = 0;
                                dst_v(i, j, k) = 2 * bnd_condition.right.y - dst_v(i - 1, j, k);
                                dst_w(i, j, k) = 2 * bnd_condition.right.z - dst_w(i - 1, j, k);
                                break;
                        }
                    }

                    if (cell_type.test(has_fluid_right))
                    {
                        switch(bnd_condition.left_type)
                        {
                            case noslip:
                                dst_u(i, j, k) = 0;
                                dst_v(i, j, k) = 2 * bnd_condition.left.y - dst_v(i + 1, j, k);
                                dst_w(i, j, k) = 2 * bnd_condition.left.z - dst_w(i + 1, j, k);
                                break;
                        }
                    }

                    if (cell_type.test(has_fluid_front))
                    {
                        switch(bnd_condition.back_type)
                        {
                            case noslip:
                                dst_u(i, j, k) = 2 * bnd_condition.back.x - dst_u(i, j - 1, k);
                                dst_v(i, j - 1, k) = 0;
                                dst_w(i, j, k) = 2 * bnd_condition.back.z - dst_w(i, j - 1, k);
                                break;
                        }
                    }

                    if (cell_type.test(has_fluid_back))
                    {
                        switch(bnd_condition.front_type)
                        {
                            case noslip:
                                dst_u(i, j, k) = 2 * bnd_condition.front.x - dst_u(i, j + 1, k);
                                dst_v(i, j, k) = 0;
                                dst_w(i, j, k) = 2 * bnd_condition.front.z - dst_w(i, j + 1, k);
                                break;
                        }
                    }
                }

            }
    };

    template<>
    struct stencils<STENCIL_COMPUTE_FG>
    {
        static void call(partition_data<Real>& dst_f, partition_data<Real>& dst_g,
            partition_data<Real>& dst_h,
            partition_data<Real> const& src_u, partition_data<Real> const& src_v,
            partition_data<Real> const& src_w,
            partition_data<std::bitset<9> > const& cell_types,
            std::vector<index> const& boundary_cells,
            std::vector<index> const& obstacle_cells,
            std::vector<index> const& fluid_cells,
            Real re, Real gx, Real gy, Real gz, Real beta, Real dx, Real dy, Real dz,
            Real dx_sq, Real dy_sq, Real dz_sq, Real dt, Real alpha
           )
        {
            for (auto const& ind : boundary_cells)
            {
                auto const i = ind.x;
                auto const j = ind.y;
                auto const k = ind.z;

                auto const& cell_type = cell_types(i, j, k);

                if (cell_type.test(has_fluid_top))
                {
                    dst_h(i, j, k) = src_w(i, j, k);
                }

                if (cell_type.test(has_fluid_right))
                {
                    dst_f(i, j, k) = src_u(i, j, k);
                }

                if (cell_type.test(has_fluid_back))
                {
                    dst_g(i, j, k) = src_v(i, j, k);
                }
            }

            for (auto const& ind : fluid_cells)
            {
                auto const i = ind.x;
                auto const j = ind.y;
                auto const k = ind.z;

                auto const& cell_type = cell_types(i, j, k);

                dst_f(i, j, k) =
                    src_u(i, j, k) + cell_type.test(has_fluid_right) *
                      dt * (
                        1. / re * (util::derivatives::second_derivative_fwd_bkwd_x(src_u, i, j, k, dx_sq)
                                    + util::derivatives::second_derivative_fwd_bkwd_y(src_u, i, j, k, dy_sq)
                                    + util::derivatives::second_derivative_fwd_bkwd_z(src_u, i, j, k, dz_sq)
                                  )

                        - util::derivatives::first_derivative_of_square_x(src_u, i, j, k, dx, alpha)
                        - util::derivatives::first_derivative_of_uv_y(src_u, src_v, i, j, k, dy, alpha)
                    //    - util::derivatives::first_derivative_of_uw_z(src_u, src_w, i, j, k, dz, alpha)
                        + gx
                    );

                dst_g(i, j, k) =
                    src_v(i, j, k) + cell_type.test(has_fluid_back) *
                    dt * (
                        1. / re * (util::derivatives::second_derivative_fwd_bkwd_x(src_v, i, j, k, dx_sq)
                                   + util::derivatives::second_derivative_fwd_bkwd_y(src_v, i, j, k, dy_sq)
                                   + util::derivatives::second_derivative_fwd_bkwd_z(src_v, i, j, k, dz_sq)
                                   )

                        - util::derivatives::first_derivative_of_uv_x(src_u, src_v, i, j, k, dx, alpha)
                        - util::derivatives::first_derivative_of_square_y(src_v, i, j, k, dy, alpha)
                      //  - util::derivatives::first_derivative_of_vw_z(src_v, src_w, i, j, k, dz, alpha)
                        + gy
                    );

                dst_h(i, j, k) =
                    src_w(i, j, k) + cell_type.test(has_fluid_top) *
                     dt * (
                        1. / re * (util::derivatives::second_derivative_fwd_bkwd_x(src_w, i, j, k, dx_sq)
                                    + util::derivatives::second_derivative_fwd_bkwd_y(src_w, i, j, k, dy_sq)
                                    + util::derivatives::second_derivative_fwd_bkwd_z(src_w, i, j, k, dz_sq)
                                  )

                        - util::derivatives::first_derivative_of_uw_x(src_u, src_w, i, j, k, dx, alpha)
                        - util::derivatives::first_derivative_of_vw_y(src_v, src_w, i, j, k, dy, alpha)
                        - util::derivatives::first_derivative_of_square_z(src_w, i, j, k, dz, alpha)
                        + gz
                    );
            }
        }
    };


    template<>
    struct stencils<STENCIL_COMPUTE_RHS>
    {
        static void call(partition_data<Real>& dst_rhs,
            partition_data<Real> const& src_f, partition_data<Real> const& src_g,
            partition_data<Real> const& src_h,
            partition_data<std::bitset<9> > const& cell_types,
            std::vector<index> const& fluid_cells,
            Real dx, Real dy, Real dz, Real dt)
        {
            for (auto const& ind : fluid_cells)
            {
                auto const i = ind.x;
                auto const j = ind.y;
                auto const k = ind.z;

                dst_rhs(i, j, k) =
                    1. / dt *   (
                                    (src_f(i, j, k) - src_f(i - 1, j, k)) / dx
                                    +
                                    (src_g(i, j, k) - src_g(i, j - 1, k)) / dy
                               //     +
                                //    (src_h(i, j, k) - src_h(i, j, k - 1)) / dz
                                );
            }

        }
    };


    template<>
    struct stencils<STENCIL_SET_P_OBSTACLE>
    {
        static void call(partition_data<Real>& dst_p,
            partition_data<std::bitset<9> > const& cell_types,
            std::vector<index> const& boundary_cells,
            std::vector<index> const& obstacle_cells,
            util::cancellation_token token)
        {
            if (!token.was_cancelled())
            {
                for (auto const& ind : boundary_cells)
                {
                    auto const i = ind.x;
                    auto const j = ind.y;
                    auto const k = ind.z;

                    auto const& cell_type = cell_types(i, j, k);

                    if (cell_type.test(has_fluid_left))
                    {
                        dst_p(i, j, k) = dst_p(i - 1, j, k);
                    }

                    if (cell_type.test(has_fluid_right))
                    {
                        dst_p(i, j, k) = dst_p(i + 1, j, k);
                    }

                    if (cell_type.test(has_fluid_bottom))
                    {
                        dst_p(i, j, k) = dst_p(i, j - 1, k);
                    }

                    if (cell_type.test(has_fluid_top))
                    {
                        dst_p(i, j, k) = dst_p(i, j + 1, k);
                    }

                    if (cell_type.test(has_fluid_front))
                    {
                        dst_p(i, j, k) = dst_p(i, j, k - 1);
                    }

                    if (cell_type.test(has_fluid_back))
                    {
                        dst_p(i, j, k) = dst_p(i, j, k + 1);
                    }
                }

            }
        }
    };

    /*
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
            {
          //       std::cout << "SOR " << iter << " | " << xd << " " << yd << std::endl;
          //  hpx::this_thread::sleep_for(boost::chrono::milliseconds(3000));

                #ifdef WITH_FOR_EACH
                hpx::parallel::for_each(
                hpx::parallel::par,
                std::begin(fluid_cells), std::end(fluid_cells),
                [&](auto const& idx_pair)
                {
                #else
                for (auto& idx_pair : fluid_cells)
                {
                #endif
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
                #ifdef WITH_FOR_EACH
                );
                #endif
            }
        }
    };
    */


    template<>
    struct stencils<STENCIL_JACOBI>
    {
        //TODO remove idx, idy (was for debug)
        static void call(partition_data<Real>& dst_p,
            partition_data<Real> const& src_rhs,
            std::vector<index> const& fluid_cells,
            Real dx_sq, Real dy_sq, Real dz_sq, util::cancellation_token token)
        {
            if (!token.was_cancelled())
            {
                for (auto const& ind : fluid_cells)
                {
                    auto const i = ind.x;
                    auto const j = ind.y;
                    auto const k = ind.z;

                    dst_p(i, j, k) =
                        (dz_sq * (dy_sq * (dst_p(i + 1, j, k) + dst_p(i - 1, j, k) - src_rhs(i, j, k) * dx_sq)
                                    + (dst_p(i, j + 1, j) + dst_p(i, j - 1, k)) * dx_sq
                                 )
                        + dx_sq * dy_sq * (dst_p(i, j, k + 1) + dst_p(i, j, k - 1))
                        )
                        /
                        (2 * (dx_sq * (dy_sq + dz_sq) + dy_sq * dz_sq));
                }
            }
        }
    };

    /*
    template<>
    struct stencils<STENCIL_TEST>
    {
        //TODO remove idx, idy (was for debug)
        static void call(partition_data<Real>& dst_p,
            partition_data<Real> const& src_rhs,
            std::vector<std::pair<std::size_t, std::size_t> > const& fluid_cells,
            std::vector<std::pair<std::size_t, std::size_t> > const& boundary_cells,
            std::vector<std::pair<std::size_t, std::size_t> > const& obstacle_cells,
            partition_data<std::bitset<6> > const& cell_types,
            Real dx_sq, Real dy_sq, util::cancellation_token token)
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
        }
    };

    */

    template<>
    struct stencils<STENCIL_COMPUTE_RESIDUAL>
    {
        static Real call(partition_data<Real> const& src_p,
            partition_data<Real> const& src_rhs,
            std::vector<index> const& fluid_cells,
            Real over_dx_sq, Real over_dy_sq, Real over_dz_sq, util::cancellation_token token)
        {
            Real local_residual = 0;

            for (auto const& ind : fluid_cells)
                {
                    auto const i = ind.x;
                    auto const j = ind.y;
                    auto const k = ind.z;

                    Real tmp =
                        (src_p(i + 1, j, k) - 2 * src_p(i, j, k) + src_p(i - 1, j, k)) / over_dx_sq
                        + (src_p(i, j + 1, k) - 2 * src_p(i, j, k) + src_p(i, j - 1, k)) / over_dy_sq
                        + (src_p(i, j, k + 1) - 2 * src_p(i, j, k) + src_p(i, j, k - 1)) / over_dz_sq
                        - src_rhs(i, j, k);

                    local_residual += std::pow(tmp, 2);
                }

            return local_residual;
        }
    };

    template<>
    struct stencils<STENCIL_UPDATE_VELOCITY>
    {
        static triple<Real> call(partition_data<Real>& dst_u,
            partition_data<Real>& dst_v,
            partition_data<Real>& dst_w,
            partition_data<Real> const& src_f,
            partition_data<Real> const& src_g,
            partition_data<Real> const& src_h,
            partition_data<Real> const& src_p,
            partition_data<std::bitset<9> > const& cell_types,
            std::vector<index> const& fluid_cells,
            Real dt, Real over_dx, Real over_dy, Real over_dz
            )
        {
            triple<Real> max_uvw(0);

            for (auto const& ind : fluid_cells)
            {
                auto const i = ind.x;
                auto const j = ind.y;
                auto const k = ind.z;

                auto const& cell_type = cell_types(i, j, k);

                if (cell_type.test(has_fluid_right))
                {
                    dst_u(i, j, k) = src_f(i, j, k) - dt * over_dx *
                        (src_p(i + 1, j, k) - src_p(i, j, k));

                    max_uvw.x = std::abs(dst_u(i, j, k)) > max_uvw.x ? std::abs(dst_u(i, j, k)) : max_uvw.x;
                }

                if (cell_type.test(has_fluid_back))
                {
                    dst_v(i, j, k) = src_g(i, j, k) - dt * over_dy *
                        (src_p(i, j + 1, k) - src_p(i, j, k));

                    max_uvw.y = std::abs(dst_v(i, j, k)) > max_uvw.y ? std::abs(dst_v(i, j, k)) : max_uvw.y;
                }

                if (cell_type.test(has_fluid_top))
                {
                 //   dst_w(i, j, k) = src_h(i, j, k) - dt * over_dz *
                 //       (src_p(i, j, k + 1) - src_p(i, j, k));

                    max_uvw.z = std::abs(dst_w(i, j, k)) > max_uvw.z ? std::abs(dst_w(i, j, k)) : max_uvw.z;
                }
            }

            return max_uvw;
        }
    };
}
}

#endif
