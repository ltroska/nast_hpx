#ifndef NAST_HPX_GRID_STENCILS_HPP_
#define NAST_HPX_GRID_STENCILS_HPP_

#include "partition_data.hpp"
#include "util/derivatives.hpp"
#include "util/cancellation_token.hpp"

#include <vector>

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
        static void call(partition_data<double>& dst, partition_data<double> const& src,
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
        static void call(partition_data<double>& dst_u, partition_data<double>& dst_v,
                partition_data<double>& dst_w,
                partition_data<std::bitset<9> > const& cell_types,
                std::vector<index>::iterator beginIt,
                std::vector<index>::iterator endIt,
                boundary_condition const& bnd_condition)
            {
                for (auto it = beginIt; it < endIt; ++it)
                {
                    auto const& ind = *it;

                    auto const i = ind.x;
                    auto const j = ind.y;
                    auto const k = ind.z;

                    auto const& cell_type = cell_types(i, j, k);

                    if (cell_type[is_boundary])
                    {
                        if (cell_type[has_fluid_bottom])
                        {
                            switch(bnd_condition.top_type)
                            {
                                case noslip:
                                    dst_u(i, j, k) = 2 * bnd_condition.top.x - dst_u(i, j, k - 1);
                                    dst_v(i, j, k) = 2 * bnd_condition.top.y - dst_v(i, j, k - 1);
                                    dst_w(i, j, k - 1) = 0;
                                    break;

                                case slip:
                                    dst_u(i, j, k) = dst_u(i, j, k - 1);
                                    dst_v(i, j, k) = dst_v(i, j, k - 1);
                                    dst_w(i, j, k - 1) = 0;
                                    break;

                                case outstream:
                                    dst_u(i, j, k) = dst_u(i, j, k - 1);
                                    dst_v(i, j, k) = dst_v(i, j, k - 1);
                                    dst_w(i, j, k - 1) = dst_w(i, j, k - 2);
                                    break;

                                case instream:
                                    dst_u(i, j, k) = 2 * bnd_condition.top.x - dst_u(i, j, k - 1);
                                    dst_v(i, j, k) = 2 * bnd_condition.top.y - dst_v(i, j, k - 1);
                                    dst_w(i, j, k - 1) = bnd_condition.top.z;
                                    break;
                            }
                        }

                        if (cell_type[has_fluid_top])
                        {
                            switch(bnd_condition.bottom_type)
                            {
                                case noslip:
                                    dst_u(i, j, k) = 2 * bnd_condition.bottom.x - dst_u(i, j, k + 1);
                                    dst_v(i, j, k) = 2 * bnd_condition.bottom.y - dst_v(i, j, k + 1);
                                    dst_w(i, j, k) = 0;
                                    break;

                                case slip:
                                    dst_u(i, j, k) = dst_u(i, j, k + 1);
                                    dst_v(i, j, k) = dst_v(i, j, k + 1);
                                    dst_w(i, j, k) = 0;
                                    break;

                                case outstream:
                                    dst_u(i, j, k) = dst_u(i, j, k + 1);
                                    dst_v(i, j, k) = dst_v(i, j, k + 1);
                                    dst_w(i, j, k) = dst_w(i, j, k + 1);
                                    break;

                                case instream:
                                    dst_u(i, j, k) = 2 * bnd_condition.bottom.x - dst_u(i, j, k + 1);
                                    dst_v(i, j, k) = 2 * bnd_condition.bottom.y - dst_v(i, j, k + 1);
                                    dst_w(i, j, k) = bnd_condition.bottom.z;
                                    break;
                            }
                        }

                        if (cell_type[has_fluid_left])
                        {
                            switch(bnd_condition.right_type)
                            {
                                case noslip:
                                    dst_u(i - 1, j, k) = 0;
                                    dst_v(i, j, k) = 2 * bnd_condition.right.y - dst_v(i - 1, j, k);
                                    dst_w(i, j, k) = 2 * bnd_condition.right.z - dst_w(i - 1, j, k);
                                    break;

                                case slip:
                                    dst_u(i - 1, j, k) = 0;
                                    dst_v(i, j, k) = dst_v(i - 1, j, k);
                                    dst_w(i, j, k) = dst_w(i - 1, j, k);
                                    break;

                                case outstream:
                                    dst_u(i - 1, j, k) = dst_u(i - 2, j, k);
                                    dst_v(i, j, k) = dst_u(i - 1, j, k);
                                    dst_w(i, j, k) = dst_u(i - 1, j, k);
                                    break;

                                case instream:
                                    dst_u(i - 1, j, k) = bnd_condition.right.x;
                                    dst_v(i, j, k) = 2 * bnd_condition.right.y - dst_v(i - 1, j, k);
                                    dst_w(i, j, k) = 2 * bnd_condition.right.z - dst_w(i - 1, j, k);
                                    break;
                            }
                        }

                        if (cell_type[has_fluid_right])
                        {
                            switch(bnd_condition.left_type)
                            {
                                case noslip:
                                    dst_u(i, j, k) = 0;
                                    dst_v(i, j, k) = 2 * bnd_condition.left.y - dst_v(i + 1, j, k);
                                    dst_w(i, j, k) = 2 * bnd_condition.left.z - dst_w(i + 1, j, k);
                                    break;

                                case slip:
                                    dst_u(i, j, k) = 0;
                                    dst_v(i, j, k) = dst_v(i + 1, j, k);
                                    dst_w(i, j, k) = dst_w(i + 1, j, k);
                                    break;

                                case outstream:
                                    dst_u(i, j, k) = dst_u(i + 1, j, k);
                                    dst_v(i, j, k) = dst_v(i + 1, j, k);
                                    dst_w(i, j, k) = dst_w(i + 1, j, k);
                                    break;

                                case instream:
                                    dst_u(i, j, k) = bnd_condition.left.x;
                                    dst_v(i, j, k) = 2 * bnd_condition.left.y - dst_v(i + 1, j, k);
                                    dst_w(i, j, k) = 2 * bnd_condition.left.z - dst_w(i + 1, j, k);
                                    break;
                            }
                        }

                        if (cell_type[has_fluid_front])
                        {
                            switch(bnd_condition.back_type)
                            {
                                case noslip:
                                    dst_u(i, j, k) = 2 * bnd_condition.back.x - dst_u(i, j - 1, k);
                                    dst_v(i, j - 1, k) = 0;
                                    dst_w(i, j, k) = 2 * bnd_condition.back.z - dst_w(i, j - 1, k);
                                    break;

                                case slip:
                                    dst_u(i, j, k) = dst_u(i, j - 1, k);
                                    dst_v(i, j - 1, k) = 0;
                                    dst_w(i, j, k) = dst_w(i, j - 1, k);
                                    break;

                                case outstream:
                                    dst_u(i, j, k) = dst_u(i, j - 1, k);
                                    dst_v(i, j - 1, k) = dst_v(i, j - 2, k);
                                    dst_w(i, j, k) = dst_w(i, j - 1, k);
                                    break;

                                case instream:
                                    dst_u(i, j, k) = 2 * bnd_condition.back.x - dst_u(i, j - 1, k);
                                    dst_v(i, j - 1, k) = bnd_condition.back.y;
                                    dst_w(i, j, k) = 2 * bnd_condition.back.z - dst_w(i, j - 1, k);
                                    break;
                            }
                        }

                        if (cell_type[has_fluid_back])
                        {
                            switch(bnd_condition.front_type)
                            {
                                case noslip:
                                    dst_u(i, j, k) = 2 * bnd_condition.front.x - dst_u(i, j + 1, k);
                                    dst_v(i, j, k) = 0;
                                    dst_w(i, j, k) = 2 * bnd_condition.front.z - dst_w(i, j + 1, k);
                                    break;

                                case slip:
                                    dst_u(i, j, k) = dst_u(i, j + 1, k);
                                    dst_v(i, j, k) = 0;
                                    dst_w(i, j, k) = dst_w(i, j + 1, k);
                                    break;

                                case outstream:
                                    dst_u(i, j, k) = dst_u(i, j + 1, k);
                                    dst_v(i, j, k) = dst_v(i, j + 1, k);
                                    dst_w(i, j, k) = dst_w(i, j + 1, k);
                                    break;

                                case instream:
                                    dst_u(i, j, k) = 2 * bnd_condition.front.x - dst_u(i, j + 1, k);
                                    dst_v(i, j, k) = bnd_condition.front.y;
                                    dst_w(i, j, k) = 2 * bnd_condition.front.z - dst_w(i, j + 1, k);
                                    break;
                            }
                        }
                    }
                    else
                    {
                        if (!cell_type[has_fluid_right] && !cell_type[has_fluid_left])
                            dst_u(i, j, k) =
                                (
                                - dst_u(i, j, k + 1) * cell_type[has_fluid_top]
                                - dst_u(i, j, k - 1) * cell_type[has_fluid_bottom]
                                - dst_u(i, j + 1, k) * cell_type[has_fluid_back]
                                - dst_u(i, j - 1, k) * cell_type[has_fluid_front]
                                )
                                /
                                (cell_type[has_fluid_top] + cell_type[has_fluid_bottom]
                                + cell_type[has_fluid_front] + cell_type[has_fluid_back]
                                );

                        if (!cell_type[has_fluid_back] && !cell_type[has_fluid_front])
                            dst_v(i, j, k) =
                                (
                                - dst_v(i, j, k + 1) * cell_type[has_fluid_top]
                                - dst_v(i, j, k - 1) * cell_type[has_fluid_bottom]
                                - dst_v(i + 1, j, k) * cell_type[has_fluid_right]
                                - dst_v(i - 1, j, k) * cell_type[has_fluid_left]
                                )
                                /
                                (cell_type[has_fluid_top] + cell_type[has_fluid_bottom]
                                + cell_type[has_fluid_right] + cell_type[has_fluid_left]
                                );

                        if (!cell_type[has_fluid_top] && !cell_type[has_fluid_bottom])
                            dst_w(i, j, k) =
                                (
                                - dst_w(i, j + 1, k) * cell_type[has_fluid_back]
                                - dst_w(i, j - 1, k) * cell_type[has_fluid_front]
                                - dst_w(i + 1, j, k) * cell_type[has_fluid_right]
                                - dst_w(i - 1, j, k) * cell_type[has_fluid_left]
                                )
                                /
                                (cell_type[has_fluid_back] + cell_type[has_fluid_front]
                                + cell_type[has_fluid_right] + cell_type[has_fluid_left]
                                );
                    }
                }

            }
    };

    template<>
    struct stencils<STENCIL_COMPUTE_FG>
    {
        static void call(partition_data<double>& dst_f, partition_data<double>& dst_g,
            partition_data<double>& dst_h,
            partition_data<double> const& src_u, partition_data<double> const& src_v,
            partition_data<double> const& src_w,
            partition_data<std::bitset<9> > const& cell_types,
            std::vector<index>::iterator beginObstacle,
            std::vector<index>::iterator endObstacle,
            std::vector<index>::iterator beginFluid,
            std::vector<index>::iterator endFluid,
            double re, double gx, double gy, double gz, double dx, double dy, double dz,
            double dx_sq, double dy_sq, double dz_sq, double dt, double alpha
           )
        {
            for (auto it = beginObstacle; it < endObstacle; ++it)
            {
                auto const& ind = *it;

                auto const i = ind.x;
                auto const j = ind.y;
                auto const k = ind.z;

                auto const& cell_type = cell_types(i, j, k);

                if (cell_type[has_fluid_top])
                {
                    dst_h(i, j, k) = src_w(i, j, k);
                }

                if (cell_type[has_fluid_right])
                {
                    dst_f(i, j, k) = src_u(i, j, k);
                }

                if (cell_type[has_fluid_back])
                {
                    dst_g(i, j, k) = src_v(i, j, k);
                }
            }

            for (auto it = beginFluid; it < endFluid; ++it)
            {
                auto const& ind = *it;

                auto const i = ind.x;
                auto const j = ind.y;
                auto const k = ind.z;

                auto const& cell_type = cell_types(i, j, k);

                dst_f(i, j, k) = src_u(i, j, k);

                if (cell_type[has_fluid_right])
                    dst_f(i, j, k) +=
                          dt * (
                            1. / re * (util::derivatives::second_derivative_fwd_bkwd_x(src_u, i, j, k, dx_sq)
                                        + util::derivatives::second_derivative_fwd_bkwd_y(src_u, i, j, k, dy_sq)
                                        + util::derivatives::second_derivative_fwd_bkwd_z(src_u, i, j, k, dz_sq)
                                      )

                            - util::derivatives::first_derivative_of_square_x(src_u, i, j, k, dx, alpha)
                            - util::derivatives::first_derivative_of_uv_y(src_u, src_v, i, j, k, dy, alpha)
                            - util::derivatives::first_derivative_of_uw_z(src_u, src_w, i, j, k, dz, alpha)
                            + gx
                        );

                dst_g(i, j, k) = src_v(i, j, k);

                if (cell_type[has_fluid_back])
                    dst_g(i, j, k) +=
                        dt * (
                            1. / re * (util::derivatives::second_derivative_fwd_bkwd_x(src_v, i, j, k, dx_sq)
                                       + util::derivatives::second_derivative_fwd_bkwd_y(src_v, i, j, k, dy_sq)
                                       + util::derivatives::second_derivative_fwd_bkwd_z(src_v, i, j, k, dz_sq)
                                       )

                            - util::derivatives::first_derivative_of_uv_x(src_u, src_v, i, j, k, dx, alpha)
                            - util::derivatives::first_derivative_of_square_y(src_v, i, j, k, dy, alpha)
                            - util::derivatives::first_derivative_of_vw_z(src_v, src_w, i, j, k, dz, alpha)
                            + gy
                        );

                dst_h(i, j, k) = src_w(i, j, k);

                if (cell_type[has_fluid_top])
                    dst_h(i, j, k) +=
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
        static void call(partition_data<double>& dst_rhs,
            partition_data<double> const& src_f, partition_data<double> const& src_g,
            partition_data<double> const& src_h,
            std::vector<index>::iterator beginIt,
            std::vector<index>::iterator endIt,
            double dx, double dy, double dz, double dt)
        {
            for (auto it = beginIt; it < endIt; ++it)
            {
                auto const& ind = *it;
                auto const i = ind.x;
                auto const j = ind.y;
                auto const k = ind.z;

                dst_rhs(i, j, k) =
                    1. / dt *   (
                                    (src_f(i, j, k) - src_f(i - 1, j, k)) / dx
                                    +
                                    (src_g(i, j, k) - src_g(i, j - 1, k)) / dy
                                    +
                                    (src_h(i, j, k) - src_h(i, j, k - 1)) / dz
                                );
            }

        }
    };


    template<>
    struct stencils<STENCIL_SET_P_OBSTACLE>
    {
        static void call(partition_data<double>& dst_p,
            partition_data<std::bitset<9> > const& cell_types,
            std::vector<index>::iterator beginIt,
            std::vector<index>::iterator endIt,
            util::cancellation_token token)
        {
            if (!token.was_cancelled())
            {
                for (auto it = beginIt; it < endIt; ++it)
                {
                    auto const& ind = *it;
                    auto const i = ind.x;
                    auto const j = ind.y;
                    auto const k = ind.z;

                    auto const& cell_type = cell_types(i, j, k);

                    dst_p(i, j, k) = (
                                            dst_p(i - 1, j, k) * cell_type[has_fluid_left]
                                          + dst_p(i + 1, j, k) * cell_type[has_fluid_right]
                                          + dst_p(i, j - 1, k) * cell_type[has_fluid_front]
                                          + dst_p(i, j + 1, k) * cell_type[has_fluid_back]
                                          + dst_p(i, j, k - 1) * cell_type[has_fluid_bottom]
                                          + dst_p(i, j, k + 1) * cell_type[has_fluid_top]
                                      )
                                      /
                                      (
                                            cell_type[has_fluid_left] + cell_type[has_fluid_right]
                                          + cell_type[has_fluid_bottom] + cell_type[has_fluid_top]
                                          + cell_type[has_fluid_front] + cell_type[has_fluid_back]
                                       );
                }
            }
        }
    };

    template<>
    struct stencils<STENCIL_SOR>
    {
        static void call(partition_data<double>& dst_p,
            partition_data<double> const& src_rhs,
            std::vector<index> const& fluid_cells,
            double part1, double part2, double dx_sq, double dy_sq, double dz_sq,
            util::cancellation_token token)
        {

            if (!token.was_cancelled())
            {
          //       std::cout << "SOR " << iter << " | " << xd << " " << yd << std::endl;
          //  hpx::this_thread::sleep_for(boost::chrono::milliseconds(3000));

                for (auto const& idx_pair : fluid_cells)
                {
                    auto i = idx_pair.x;
                    auto j = idx_pair.y;
                    auto k = idx_pair.z;

                    dst_p(i, j, k) =
                        part1 * dst_p(i, j, k)
                        + part2 * (
                                (dst_p(i + 1, j, k) + dst_p(i - 1, j, k)) / dx_sq
                                + (dst_p(i, j + 1, k) + dst_p(i, j - 1, k)) / dy_sq
                                + (dst_p(i, j, k + 1) + dst_p(i, j, k - 1)) / dz_sq
                                - src_rhs(i, j, k)
                        );
                }
            }
        }
    };

    template<>
    struct stencils<STENCIL_JACOBI>
    {
        static void call(partition_data<double>& dst_p,
            partition_data<double> const& src_rhs,
            std::vector<index>::iterator beginIt,
            std::vector<index>::iterator endIt,
            double dx_sq, double dy_sq, double dz_sq, util::cancellation_token token)
        {
            if (!token.was_cancelled())
            {
                for (auto it = beginIt; it < endIt; ++it)
                {
                    auto const& ind = *it;
                    auto const i = ind.x;
                    auto const j = ind.y;
                    auto const k = ind.z;

                    dst_p(i, j, k) =
                        dx_sq * dy_sq * dz_sq / (2. * (dx_sq * dy_sq + dx_sq * dz_sq + dy_sq * dz_sq)) *
                        (
                                (dst_p(i + 1, j, k) + dst_p(i - 1, j, k)) / dx_sq
                                + (dst_p(i, j + 1, k) + dst_p(i, j - 1, k)) / dy_sq
                                + (dst_p(i, j, k + 1) + dst_p(i, j, k - 1)) / dz_sq
                                - src_rhs(i, j, k)
                        );
                }
            }
        }
    };

    template<>
    struct stencils<STENCIL_COMPUTE_RESIDUAL>
    {
        static double call(partition_data<double> const& src_p,
            partition_data<double> const& src_rhs,
            std::vector<index>::iterator beginIt,
            std::vector<index>::iterator endIt,
            double over_dx_sq, double over_dy_sq, double over_dz_sq, util::cancellation_token token)
        {
            double local_residual = 0;

            if (!token.was_cancelled())
                for (auto it = beginIt; it < endIt; ++it)
                {
                    auto const& ind = *it;
                    auto const i = ind.x;
                    auto const j = ind.y;
                    auto const k = ind.z;

                    double tmp =
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
        static triple<double> call(partition_data<double>& dst_u,
            partition_data<double>& dst_v,
            partition_data<double>& dst_w,
            partition_data<double> const& src_f,
            partition_data<double> const& src_g,
            partition_data<double> const& src_h,
            partition_data<double> const& src_p,
            partition_data<std::bitset<9> > const& cell_types,
            std::vector<index>::iterator beginIt,
            std::vector<index>::iterator endIt,
            double dt, double over_dx, double over_dy, double over_dz
            )
        {
            triple<double> max_uvw(0);

            for (auto it = beginIt; it < endIt; ++it)
            {
                auto const& ind = *it;
                auto const i = ind.x;
                auto const j = ind.y;
                auto const k = ind.z;

                auto const& cell_type = cell_types(i, j, k);

                if (cell_type[has_fluid_right])
                {
                    dst_u(i, j, k) = src_f(i, j, k) - dt * over_dx *
                        (src_p(i + 1, j, k) - src_p(i, j, k));

                    max_uvw.x = std::abs(dst_u(i, j, k)) > max_uvw.x ? std::abs(dst_u(i, j, k)) : max_uvw.x;
                }

                if (cell_type[has_fluid_back])
                {
                    dst_v(i, j, k) = src_g(i, j, k) - dt * over_dy *
                        (src_p(i, j + 1, k) - src_p(i, j, k));

                    max_uvw.y = std::abs(dst_v(i, j, k)) > max_uvw.y ? std::abs(dst_v(i, j, k)) : max_uvw.y;
                }

                if (cell_type[has_fluid_top])
                {
                    dst_w(i, j, k) = src_h(i, j, k) - dt * over_dz *
                        (src_p(i, j, k + 1) - src_p(i, j, k));

                    max_uvw.z = std::abs(dst_w(i, j, k)) > max_uvw.z ? std::abs(dst_w(i, j, k)) : max_uvw.z;
                }
            }

            return max_uvw;
        }
    };
}
}

#endif
