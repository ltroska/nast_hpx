#ifndef NAST_HPX_GRID_STENCILS_HPP
#define NAST_HPX_GRID_STENCILS_HPP

#include <vector>

#include "partition_data.hpp"
#include "fd_stencils.hpp"

namespace nast_hpx { namespace grid {
    
    static const std::size_t STENCIL_NONE = 20;
    static const std::size_t STENCIL_SET_VELOCITY_BOUNDARY_NOSLIP = 21;
    static const std::size_t STENCIL_SET_VELOCITY_BOUNDARY_SLIP = 22;
    static const std::size_t STENCIL_SET_VELOCITY_BOUNDARY_OUTSTREAM = 23;
    static const std::size_t STENCIL_SET_VELOCITY_BOUNDARY_INSTREAM = 24;
    static const std::size_t STENCIL_SET_VELOCITY_OBSTACLE = 25;
    static const std::size_t STENCIL_SET_VELOCITY = 26;
    static const std::size_t STENCIL_COMPUTE_FG = 27;
    static const std::size_t STENCIL_COMPUTE_RHS = 28;
    static const std::size_t STENCIL_SET_P = 29;
    static const std::size_t STENCIL_SOR_CYCLE = 30;
    static const std::size_t STENCIL_SOR_RESIDUAL = 31;
    static const std::size_t STENCIL_UPDATE_VELOCITY = 32;
    
    typedef std::pair<std::size_t, std::size_t> range_type;
    
    namespace detail
    {
        namespace velocity
        {
            template<direction dir>
            void set_noslip(partition_data<Real>& dst_u, partition_data<Real>& dst_v,
                partition_data<Real> const& src_u, partition_data<Real> const& src_v,
                range_type indices, Real u_boundary, Real v_boundary);
            
            template<>
            void set_noslip<LEFT>(partition_data<Real>& dst_u, partition_data<Real>& dst_v,
                partition_data<Real> const& src_u, partition_data<Real> const& src_v,
                range_type indices, Real u_boundary, Real v_boundary)
            {
                constexpr std::size_t x = 1;
                for (std::size_t y = indices.first; y < indices.second; ++y)
                {
                    dst_u(x, y) = 0;
                    dst_v(x, y) = -src_v(x + 1, y);
                }
                    
            }        
            
            template<>
            void set_noslip<RIGHT>(partition_data<Real>& dst_u, partition_data<Real>& dst_v,
                partition_data<Real> const& src_u, partition_data<Real> const& src_v,
                range_type indices, Real u_boundary, Real v_boundary)
            {
                const std::size_t x = src_u.size_x_ - 2;
                for (std::size_t y = indices.first; y < indices.second; ++y)
                {
                    dst_u(x - 1, y) = 0;
                    dst_v(x, y) = -src_v(x - 1, y);
                }
                    
            }        
            
            template<>
            void set_noslip<BOTTOM>(partition_data<Real>& dst_u, partition_data<Real>& dst_v,
                partition_data<Real> const& src_u, partition_data<Real> const& src_v,
                range_type indices, Real u_boundary, Real v_boundary)
            {
                constexpr std::size_t y = 1;
                for (std::size_t x = indices.first; x < indices.second; ++x)
                {
                    dst_u(x, y) = -src_u(x, y + 1);
                    dst_v(x, y) = 0;
                }
                    
            }    
        
            template<>
            void set_noslip<TOP>(partition_data<Real>& dst_u, partition_data<Real>& dst_v,
                partition_data<Real> const& src_u, partition_data<Real> const& src_v,
                range_type indices, Real u_boundary, Real v_boundary)
            {
                const std::size_t y = src_u.size_y_ - 2;
                for (std::size_t x = indices.first; x < indices.second; ++x)
                {
                    dst_u(x, y) = 2 * u_boundary - src_u(x, y - 1);
                    dst_v(x, y - 1) = 0;
                }
                    
            }
        }
    }
    
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
    struct stencils<STENCIL_SET_VELOCITY_BOUNDARY_NOSLIP>
    {
        template<direction dir>
        static void call(partition_data<Real>& dst_u, partition_data<Real>& dst_v,
            partition_data<Real> const& src_u, partition_data<Real> const& src_v,
            range_type indices, Real u_boundary, Real v_boundary)
            {
                detail::velocity::set_noslip<dir>(dst_u, dst_v, src_u, src_v, indices, u_boundary, v_boundary);
            }                
    };      
    
    template<>
    struct stencils<STENCIL_SET_VELOCITY_OBSTACLE>
    {
        static void call(partition_data<Real>& dst_u, partition_data<Real>& dst_v,
            partition_data<Real> const& src_u, partition_data<Real> const& src_v,
            partition_data<std::bitset<6> > const& cell_types,
            range_type x_range, range_type y_range)
            {
                for (std::size_t y = y_range.first; y < y_range.second; ++y)
                    for (std::size_t x = x_range.first; x < x_range.second; ++x)
                    {
                        auto const& cell_type = cell_types(x, y);
                        
                        if (cell_type.test(is_fluid))
                        {      
                            // only set fluid cells adjacent to obstacle
                            // fluid adjacent to boundary is set in boundary
                            // stencil
                            if (!cell_type.test(has_fluid_east) && cell_types(x + 1, y).test(has_fluid_west))
                                dst_u(x, y) = 0;
                  
                            if (!cell_type.test(has_fluid_north) && cell_types(x, y + 1).test(has_fluid_south))
                                dst_v(x, y) = 0;
                        }
                        else
                        {
                            dst_u(x, y) =
                                - src_u(x, y - 1) * cell_type.test(has_fluid_south)
                                - src_u(x, y + 1) * cell_type.test(has_fluid_north);
                                
                            dst_v(x, y) =
                                - src_v(x - 1, y) * cell_type.test(has_fluid_west)
                                - src_v(x + 1, y) * cell_type.test(has_fluid_east); 
                        }
                    }                
            }                
    };

    template<>
    struct stencils<STENCIL_SET_VELOCITY>
    {
        static void call(partition_data<Real>& dst_u, partition_data<Real>& dst_v,
            partition_data<Real> const& src_u, partition_data<Real> const& src_v,
            partition_data<std::bitset<6> > const& cell_types,
            boundary_data u_bnd, boundary_data v_bnd, range_type x_range,
            range_type y_range)
            {                
                for (std::size_t y = y_range.first; y < y_range.second; ++y)
                    for (std::size_t x = x_range.first; x < x_range.second; ++x)
                    {
                        auto const& cell_type = cell_types(x, y);
                        
                        if (cell_type.test(is_boundary) /*&& !cell_type.test(is_fluid)*/)
                        {
                            if (cell_type.test(has_fluid_east))
                            {
                                dst_u(x, y) = 0;
                                dst_v(x, y) = -src_v(x + 1, y); 
                            }
                            
                            else if (cell_type.test(has_fluid_west))
                            {
                                dst_u(x - 1, y) = 0;
                                dst_v(x, y) = -src_v(x - 1, y);
                            }
                            
                            else if (cell_type.test(has_fluid_north))
                            {
                                dst_u(x, y) = -src_u(x, y + 1);
                                dst_v(x, y) = 0;
                            }
                            
                            else if (cell_type.test(has_fluid_south))
                            {
                                dst_u(x, y) = 2 * u_bnd.top - src_u(x, y - 1);
                                dst_v(x, y - 1) = 0;
                            }                            
                        }

                        else if (!cell_type.test(is_fluid) && !cell_type.none())
                        {
                            dst_u(x, y) =
                                - src_u(x, y - 1) * cell_type.test(has_fluid_south)
                                - src_u(x, y + 1) * cell_type.test(has_fluid_north);
                                
                            dst_v(x, y) =
                                - src_v(x - 1, y) * cell_type.test(has_fluid_west)
                                - src_v(x + 1, y) * cell_type.test(has_fluid_east); 
                        }
                        
                        else
                        {
                            // only set fluid cells adjacent to obstacle
                            // fluid adjacent to boundary is set above
                            if (!cell_type.test(has_fluid_east) && !cell_types(x + 1, y).test(is_boundary))
                                dst_u(x, y) = 0;
                  
                            if (!cell_type.test(has_fluid_north) && !cell_types(x, y + 1).test(is_boundary))
                                dst_v(x, y) = 0;
                        }
                    }
            }                
    }; 

    template<>
    struct stencils<STENCIL_COMPUTE_FG>
    {
        static void call(partition_data<Real>& dst_f, partition_data<Real>& dst_g,
            partition_data<Real> const& src_u, partition_data<Real> const& src_v,
            partition_data<std::bitset<6> > const& cell_types,
            Real re, Real gx, Real gy, Real beta, Real dx, Real dy, Real dt,
            Real alpha, range_type x_range, range_type y_range)
        {
            for (std::size_t y = y_range.first; y < y_range.second; ++y)
                for (std::size_t x = x_range.first; x < x_range.second; ++x)
                {
                    auto const& cell_type = cell_types(x, y);
                    
                    if (!cell_type.test(is_fluid))
                    {
                        if (cell_type.test(has_fluid_east))
                            dst_f(x, y) = src_u(x, y);

                        if (cell_type.test(has_fluid_north))
                            dst_g(x, y) = src_v(x, y);
                    }
                    else
                    {                        
                        dst_f(x, y) = src_u(x, y) + cell_type.test(has_fluid_east)
                            * ( dt * (
                                1. / re * (
                                    second_derivative_fwd_bkwd_x(
                                        src_u(x + 1, y), src_u(x, y), src_u(x - 1, y), dx
                                    )
                                    + second_derivative_fwd_bkwd_y(
                                        src_u(x, y + 1), src_u(x, y), src_u(x, y - 1), dy
                                    )
                                )
                                
                                - first_derivative_of_square_x(
                                    src_u(x + 1, y), src_u(x, y), src_u(x - 1, y),
                                    dx, alpha
                                )
                                - first_derivative_of_product_y(
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
                                    second_derivative_fwd_bkwd_x(
                                        src_v(x + 1, y), src_v(x, y), src_v(x - 1, y), dx
                                    )
                                    + second_derivative_fwd_bkwd_y(
                                        src_v(x, y + 1), src_v(x, y), src_v(x, y - 1), dy
                                    )
                                )
                                
                                - first_derivative_of_product_x(
                                    src_u(x - 1, y), src_u(x, y), src_u(x, y + 1),
                                    src_u(x - 1, y + 1), src_v(x - 1, y),
                                    src_v(x, y), src_v(x + 1, y), dx, alpha
                                )
                                - first_derivative_of_square_y(
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
        }
    };
    
    template<>
    struct stencils<STENCIL_COMPUTE_RHS>
    {
        static void call(partition_data<Real>& dst_rhs,
            partition_data<Real> const& src_f, partition_data<Real> const& src_g,
            partition_data<std::bitset<6> > const& cell_types,
            Real dx, Real dy, Real dt, range_type x_range, range_type y_range)
        {
            for (std::size_t y = y_range.first; y < y_range.second; ++y)
                for (std::size_t x = x_range.first; x < x_range.second; ++x)
                {
                    auto const& cell_type = cell_types(x, y);
                    
                    if (cell_type.test(is_fluid))
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
            std::size_t idx, std::size_t idy,
            range_type x_range, range_type y_range)
        {
         //   std::cout << "set p " << idx << " " << idy << std::endl;
            
            for (std::size_t y = y_range.first; y < y_range.second; ++y)
                for (std::size_t x = x_range.first; x < x_range.second; ++x)
                {
                    auto const& cell_type = cell_types(x, y);
                    
                    if (!cell_type.test(is_fluid) && !cell_type.none())
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
    };    
        
    template<>
    struct stencils<STENCIL_SOR_CYCLE>
    {
        //TODO remove idx, idy (was for debug)
        static void call(partition_data<Real>& dst_p,
            partition_data<Real> const& src_rhs,
            partition_data<std::bitset<6> > const& cell_types,
            Real part1, Real part2, Real dx_sq, Real dy_sq,
            std::size_t idx, std::size_t idy, std::size_t iter,
            range_type x_range, range_type y_range)
        {
           /* std::cout << "sor in iteration " << iter << " " << idx << " " << idy << std::endl;
            if (hpx::get_locality_id() == 0)
                hpx::this_thread::sleep_for(boost::chrono::milliseconds(10000));
            else
            hpx::this_thread::sleep_for(boost::chrono::milliseconds(10000));*/

            for (std::size_t y = y_range.first; y < y_range.second; ++y)
                for (std::size_t x = x_range.first; x < x_range.second; ++x)
                {
                    auto const& cell_type = cell_types(x, y);
                    
                    if (cell_type.test(is_fluid))
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
    struct stencils<STENCIL_SOR_RESIDUAL>
    {
        static Real call(partition_data<Real> const& src_p,
            partition_data<Real> const& src_rhs,
            partition_data<std::bitset<6> > const& cell_types,
            Real over_dx_sq, Real over_dy_sq,
            range_type x_range, range_type y_range)
        {
           // std::cout << "sor in " << idx << " " << idy << std::endl;
            

            Real local_residual = 0;

            for (std::size_t y = y_range.first; y < y_range.second; ++y)
                for (std::size_t x = x_range.first; x < x_range.second; ++x)
                {
                    auto const& cell_type = cell_types(x, y);
                    
                    if (cell_type.test(is_fluid))
                    {
                        Real tmp =
                            (src_p(x + 1, y) - 2 * src_p(x, y) + src_p(x - 1, y)) * over_dx_sq
                            + (src_p(x, y + 1) - 2 * src_p(x, y) + src_p(x, y - 1)) * over_dy_sq
                            - src_rhs(x, y);
                                                
                        local_residual += std::pow(tmp, 2);
                    }
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
            Real dt, Real over_dx, Real over_dy,
            range_type x_range, range_type y_range)
        {
           // std::cout << "sor in " << idx << " " << idy << std::endl;
            
            Real max_u = 0;
            Real max_v = 0;
            
            for (std::size_t y = y_range.first; y < y_range.second; ++y)
                for (std::size_t x = x_range.first; x < x_range.second; ++x)
                {
                    auto const& cell_type = cell_types(x, y);
                    
                    if (cell_type.test(is_fluid))
                    {
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
                                
                            max_v = std::abs(dst_v(x, y)) > max_v ? std::abs(dst_u(x, y)) : max_v;
                        }
                    }
                }

            return std::make_pair(max_u, max_v);
        }
    };
}
}

#endif
