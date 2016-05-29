#ifndef NAST_HPX_GRID_STENCILS_HPP
#define NAST_HPX_GRID_STENCILS_HPP

#include <vector>

#include "partition_data.hpp"

namespace nast_hpx { namespace grid {
    
    static const std::size_t STENCIL_NONE = 20;
    static const std::size_t STENCIL_SET_VELOCITY_BOUNDARY_NOSLIP = 21;
    static const std::size_t STENCIL_SET_VELOCITY_BOUNDARY_SLIP = 22;
    static const std::size_t STENCIL_SET_VELOCITY_BOUNDARY_OUTSTREAM = 23;
    static const std::size_t STENCIL_SET_VELOCITY_BOUNDARY_INSTREAM = 24;
    static const std::size_t STENCIL_SET_VELOCITY_OBSTACLE = 25;
    static const std::size_t STENCIL_SET_VELOCITY = 26;
    static const std::size_t STENCIL_COMPUTE_FG = 27;
    
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

                        else if (!cell_type.test(is_fluid))
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
                            if (!cell_type.test(has_fluid_east) && cell_types(x + 1, y).test(has_fluid_west))
                                dst_u(x, y) = 0;
                  
                            if (!cell_type.test(has_fluid_north) && cell_types(x, y + 1).test(has_fluid_south))
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
                    
                    
                    
                  /*  middle_fg.first =
                        middle_uv.first
                        + type.test(3) * (dt * (
                                            1./re
                                            *   (second_derivative_fwd_bkwd_x(right_uv.first, 
                                                    middle_uv.first, left_uv.first, dx)
                                                + second_derivative_fwd_bkwd_y(top_uv.first,
                                                    middle_uv.first, bottom_uv.first, dy))

                                                - first_derivative_of_square_x(right_uv.first,
                                                    middle_uv.first, left_uv.first, dx, alpha)
                                                - first_derivative_of_product_y(right_uv.second,
                                                    middle_uv.second, bottom_uv.second,
                                                    bottomright_uv.second, bottom_uv.first,
                                                    middle_uv.first, top_uv.first, dy, alpha)

                                                + gx
                                            )
                                           /* - beta * dt / 2.
                                              * (middle_temperature 
                                                    + right_temperature )
                                              * gx
                                        );
                    middle_fg.second =
                        middle_uv.second
                            + type.test(0) * (dt * (
                                            1./re
                                            *   (second_derivative_fwd_bkwd_x(right_uv.second,
                                                    middle_uv.second, left_uv.second, dx)
                                                + second_derivative_fwd_bkwd_y(top_uv.second,
                                                    middle_uv.second, bottom_uv.second, dy))

                                            - first_derivative_of_product_x(left_uv.first,
                                                middle_uv.first, top_uv.first, topleft_uv.first,
                                                left_uv.second, middle_uv.second,
                                                right_uv.second, dx, alpha)
                                
                                            - first_derivative_of_square_y(top_uv.second,
                                                middle_uv.second, bottom_uv.second, dy, alpha)

                                            + gy
                                            )
                                            - beta * dt / 2.
                                                * (middle_temperature 
                                                    + top_temperature )
                                                * gy
                                            );*/
                }     
            
        }
    };
}
}

#endif
