#ifndef NAST_HPX_GRID_STENCILS_HPP
#define NAST_HPX_GRID_STENCILS_HPP

#include <vector>

#ifdef WITH_FOR_EACH
#include <hpx/parallel/algorithms/for_each.hpp>
#include <hpx/parallel/algorithms/reduce.hpp>
#include <hpx/parallel/algorithms/transform_reduce.hpp>
#endif

#include "partition_data.hpp"
#include "util/fd_stencils.hpp"
#include "util/cancellation_token.hpp"

namespace nast_hpx { namespace grid {

    static const std::size_t STENCIL_NONE = 20;
    static const std::size_t STENCIL_SET_VELOCITY_BOUNDARY = 21;
    static const std::size_t STENCIL_SET_VELOCITY_OBSTACLE = 22;
    static const std::size_t STENCIL_SET_VELOCITY = 23;
    static const std::size_t STENCIL_COMPUTE_FG_BOUNDARY_AND_OBSTACLE = 24;
    static const std::size_t STENCIL_COMPUTE_FG_FLUID = 25;
    static const std::size_t STENCIL_COMPUTE_FG = 33;
    static const std::size_t STENCIL_COMPUTE_RHS = 26;
    static const std::size_t STENCIL_SET_P = 27;
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
            boundary_data u_bnd, boundary_data v_bnd, boundary_type bnd_type
            )
            {
                #ifdef WITH_FOR_EACH
                hpx::parallel::for_each(
                hpx::parallel::par,
                std::begin(boundary_cells), std::end(boundary_cells),
                [&](auto const& idx_pair)
                {
                #else
                for (auto const& idx_pair : boundary_cells)
                {
                #endif
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
                                dst_v(x, y) = 2 * v_bnd.left - src_v(x + 1, y);
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
                #ifdef WITH_FOR_EACH
                );
                #endif
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
                #ifdef WITH_FOR_EACH
                hpx::parallel::for_each(
                hpx::parallel::par,
                std::begin(obstacle_cells), std::end(obstacle_cells),
                [&](auto const& idx_pair)
                {
                #else
                for (auto const& idx_pair : obstacle_cells)
                {
                #endif
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
                #ifdef WITH_FOR_EACH
                );
                #endif
            }
    };

    template<>
    struct stencils<STENCIL_SET_VELOCITY>
    {
        static void call(partition_data<Real>& dst_u, partition_data<Real>& dst_v,
            partition_data<Real> const& src_u, partition_data<Real> const& src_v,
            partition_data<std::bitset<6> > const& cell_types,
            std::vector<std::pair<std::size_t, std::size_t> > const& boundary_cells,
            std::vector<std::pair<std::size_t, std::size_t> > const& obstacle_cells,
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
                                dst_v(x, y) = 2 * v_bnd.left - src_v(x + 1, y);
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
    struct stencils<STENCIL_COMPUTE_FG_BOUNDARY_AND_OBSTACLE>
    {
        static void call(partition_data<Real>& dst_f, partition_data<Real>& dst_g,
            partition_data<Real> const& src_u, partition_data<Real> const& src_v,
            partition_data<std::bitset<6> > const& cell_types,
            std::vector<std::pair<std::size_t, std::size_t> > const& boundary_cells,
            std::vector<std::pair<std::size_t, std::size_t> > const& obstacle_cells
           )
        {
            #ifdef WITH_FOR_EACH
            auto f1 = hpx::parallel::for_each(
            hpx::parallel::par(hpx::parallel::task),
            std::begin(boundary_cells), std::end(boundary_cells),
            [&](auto const& idx_pair)
            {
            #else
            for (auto const& idx_pair : boundary_cells)
            {
            #endif
                auto const x = idx_pair.first;
                auto const y = idx_pair.second;

                auto const& cell_type = cell_types(x, y);

                if (cell_type.test(has_fluid_east))
                    dst_f(x, y) = src_u(x, y);

                if (cell_type.test(has_fluid_north))
                    dst_g(x, y) = src_v(x, y);
            }
            #ifdef WITH_FOR_EACH
            );

            auto f2 = hpx::parallel::for_each(
            hpx::parallel::par(hpx::parallel::task),
            std::begin(obstacle_cells), std::end(obstacle_cells),
            [&](auto const& idx_pair)
            {
            #else
            for (auto const& idx_pair : obstacle_cells)
            {
            #endif
                auto const x = idx_pair.first;
                auto const y = idx_pair.second;

                auto const& cell_type = cell_types(x, y);

                if (cell_type.test(has_fluid_east))
                    dst_f(x, y) = src_u(x, y);

                if (cell_type.test(has_fluid_north))
                    dst_g(x, y) = src_v(x, y);
            }
            #ifdef WITH_FOR_EACH
            );

            hpx::wait_all(f1, f2);
            #endif
        }
    };

    template<>
    struct stencils<STENCIL_COMPUTE_FG>
    {
        static void call(partition_data<Real>& dst_f, partition_data<Real>& dst_g,
            partition_data<Real> const& src_u, partition_data<Real> const& src_v,
            partition_data<std::bitset<6> > const& cell_types,
            std::vector<std::pair<std::size_t, std::size_t> > const& boundary_cells,
            std::vector<std::pair<std::size_t, std::size_t> > const& obstacle_cells,
            std::vector<std::pair<std::size_t, std::size_t> > const& fluid_cells,
            Real re, Real gx, Real gy, Real beta, Real dx, Real dy, Real dt,
            Real alpha
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
    struct stencils<STENCIL_COMPUTE_FG_FLUID>
    {
        static void call(partition_data<Real>& dst_f, partition_data<Real>& dst_g,
            partition_data<Real> const& src_u, partition_data<Real> const& src_v,
            partition_data<std::bitset<6> > const& cell_types,
            std::vector<std::pair<std::size_t, std::size_t> > const& fluid_cells,
            Real re, Real gx, Real gy, Real beta, Real dx, Real dy, Real dt,
            Real alpha)
        {
            #ifdef WITH_FOR_EACH
            hpx::parallel::for_each(
            hpx::parallel::par,
            std::begin(fluid_cells), std::end(fluid_cells),
            [&](auto const& idx_pair)
            {
            #else
            for (auto const& idx_pair : fluid_cells)
            {
            #endif
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
            #ifdef WITH_FOR_EACH
            );
            #endif
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
            #ifdef WITH_FOR_EACH
            hpx::parallel::for_each(
            hpx::parallel::par,
            std::begin(fluid_cells), std::end(fluid_cells),
            [&](auto const& idx_pair)
            {
            #else
            for (auto const& idx_pair : fluid_cells)
            {
            #endif
                auto const x = idx_pair.first;
                auto const y = idx_pair.second;

                dst_rhs(x, y) =
                    1. / dt *   (
                                    (src_f(x, y) - src_f(x - 1, y)) / dx
                                    +
                                    (src_g(x, y) - src_g(x, y - 1)) / dy
                                );
            }
            #ifdef WITH_FOR_EACH
            );
            #endif
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
                #ifdef WITH_FOR_EACH
                auto f1 = hpx::parallel::for_each(
                hpx::parallel::par(hpx::parallel::task),
                std::begin(boundary_cells), std::end(boundary_cells),
                [&](auto const& idx_pair)
                {
                #else
                for (auto& idx_pair : boundary_cells)
                {
                #endif
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
                #ifdef WITH_FOR_EACH
                );

                auto f2 = hpx::parallel::for_each(
                hpx::parallel::par(hpx::parallel::task),
                std::begin(obstacle_cells), std::end(obstacle_cells),
                [&](auto const& idx_pair)
                {
                #else
                for (auto& idx_pair : obstacle_cells)
                {
                #endif
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
                #ifdef WITH_FOR_EACH
                );

                hpx::wait_all(f1, f2);
                #endif
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
            {
               // std::cout << "jacobi " << iter << " | " << xd << " " << yd << std::endl;
              //  hpx::this_thread::sleep_for(boost::chrono::milliseconds(3000));

              //  if (hpx::get_locality_id() == 1)
              //      hpx::this_thread::sleep_for(boost::chrono::milliseconds(4000));
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
                         ( (dst_p(x + 1, y) + dst_p(x - 1, y)) * dy_sq
                            + (dst_p(x, y + 1) + dst_p(x, y - 1)) * dx_sq
                            - dx_sq * dy_sq * src_rhs(x, y))
                        /
                        (2 * (dx_sq + dy_sq));
                }
                #ifdef WITH_FOR_EACH
                );
                #endif
            }
        }
    };

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
            #ifdef WITH_FOR_EACH
                local_residual =
                hpx::parallel::transform_reduce(
                hpx::parallel::par,
                std::begin(fluid_cells), std::end(fluid_cells),
                [&](auto const& idx_pair) -> Real
                {
                    auto x = idx_pair.first;
                    auto y = idx_pair.second;

                    Real tmp =
                            (src_p(x + 1, y) - 2 * src_p(x, y) + src_p(x - 1, y)) * over_dx_sq
                            + (src_p(x, y + 1) - 2 * src_p(x, y) + src_p(x, y - 1)) * over_dy_sq
                            - src_rhs(x, y);

                    return std::pow(tmp, 2);
                }
                , 0.
                ,[](Real const a, Real const b) -> Real
                {
                    return a + b;
                }
                );

            #else
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
            #endif

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
            #ifdef WITH_FOR_EACH
            std::pair<Real, Real> max_uv =
                hpx::parallel::transform_reduce(
                    hpx::parallel::par,
                    std::begin(fluid_cells), std::end(fluid_cells),
                    [&](auto const& idx_pair) -> std::pair<Real, Real>
                    {
                        auto const x = idx_pair.first;
                        auto const y = idx_pair.second;

                        auto const& cell_type = cell_types(x, y);


                        if (cell_type.test(has_fluid_east))
                            dst_u(x, y) = src_f(x, y) - dt * over_dx *
                                (src_p(x + 1, y) - src_p(x, y));

                        if (cell_type.test(has_fluid_north))
                            dst_v(x, y) = src_g(x, y) - dt * over_dy *
                                (src_p(x, y + 1) - src_p(x, y));

                        return std::make_pair(std::abs(dst_u(x, y)), std::abs(dst_v(x, y)));
                    },
                    std::make_pair(0., 0.),
                    [](auto const& a, auto const& b) -> std::pair<Real, Real>
                    {
                        return std::make_pair(a.first > b.first ? a.first : b.first, a.second > b.second ? a.second : b.second);
                    }
                );

            return max_uv;

            #else
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
            #endif
        }
    };
}
}

#endif
