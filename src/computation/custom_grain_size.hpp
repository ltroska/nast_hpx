/** Strategy that implements the algorithm on the grid using a division
 *  into smaller sub blocks of specified size, with heavy use of dataflow
 *  to chain the computation together.
 */

#ifndef CUSTOM_GRAIN_SIZE_HPP
#define CUSTOM_GRAIN_SIZE_HPP

#include <hpx/hpx.hpp>

#include "grid/types.hpp"
#include "util/boundary_data.hpp"

namespace computation {

class custom_grain_size {

public:
    static vector_partition set_velocity_for_boundary_and_obstacles(
            vector_partition const& middle_partition,
            vector_partition const& left_partition,
            vector_partition const& right_partition,
            vector_partition const& bottom_partition,
            vector_partition const& top_partition,
            std::vector<std::vector<std::pair<uint, uint> > > const& boundary,
            std::vector<std::pair<uint, uint> > const& obstacle,
            std::vector<std::pair<uint, uint> > const& fluid,
            std::vector<std::bitset<5> > const& flag_data,
            boundary_data const& boundary_data_type,
            boundary_data const& u_boundary_data,
            boundary_data const& v_boundary_data);
    
    static scalar_partition set_temperature_for_boundary_and_obstacles(
            scalar_partition const& middle_partition,
            scalar_partition const& left_partition,
            scalar_partition const& right_partition,
            scalar_partition const& bottom_partition,
            scalar_partition const& top_partition,
            std::vector<std::vector<std::pair<uint, uint> > > const& boundary,
            boundary_data const& boundary_data_type,
            boundary_data const& temperature_boundary_data,
            uint global_i, uint global_j,
            RealType dx, RealType dy);
    
    static vector_partition compute_fg_on_fluid_cells(
        vector_partition const& middle_fg,
        vector_partition const& middle_uv,
        vector_partition const& left_uv, vector_partition const& right_uv,
        vector_partition const& bottom_uv, vector_partition const& top_uv,
        vector_partition const& bottomright_uv,
        vector_partition const& topleft_uv,
        scalar_partition const& middle_temperature,
        scalar_partition const& right_temperature,
        scalar_partition const& top_temperature,
        std::vector<std::vector<std::pair<uint, uint> > > const& boundary,
        std::vector<std::pair<uint, uint> > const& obstacle,
        std::vector<std::pair<uint, uint> > const& fluid,
        std::vector<std::bitset<5> > const& flag_data,
        RealType re, RealType gx, RealType gy, RealType beta,
        RealType dx, RealType dy, RealType dt, RealType alpha);
    
    static scalar_partition compute_temperature_on_fluid_cells(
        scalar_partition const& middle_temperature,
        scalar_partition const& left_temperature,
        scalar_partition const& right_temperature,
        scalar_partition const& bottom_temperature,
        scalar_partition const& top_temperature,
        vector_partition const& middle_uv,
        vector_partition const& left_uv,
        vector_partition const& bottom_uv,
        std::vector<std::vector<std::pair<uint, uint> > > const& boundary,
        std::vector<std::pair<uint, uint> > const& obstacle,
        std::vector<std::pair<uint, uint> > const& fluid,
        RealType re, RealType pr, RealType dx, RealType dy, RealType dt,
        RealType alpha);
    
    static scalar_partition compute_right_hand_side_on_fluid_cells(
        scalar_partition const& middle_rhs, vector_partition const& middle_fg,
        vector_partition const& left_fg, vector_partition const& bottom_fg,
        std::vector<std::pair<uint, uint> > const& fluid,
        RealType dx, RealType dy, RealType dt);
    
    static scalar_data set_pressure_for_boundary_and_obstacles(
       hpx::shared_future<scalar_data> middle_p,
        hpx::shared_future<scalar_data> left_p, hpx::shared_future<scalar_data> right_p,
        hpx::shared_future<scalar_data> bottom_p, hpx::shared_future<scalar_data> top_p,
        std::vector<std::vector<std::pair<uint, uint> > > const& boundary,
        std::vector<std::pair<uint, uint> > const& obstacle,
        std::vector<std::bitset<5> > const& flag_data);
    
    static scalar_data sor_cycle(hpx::shared_future<scalar_data> middle_p,
        hpx::shared_future<scalar_data> left_p, hpx::shared_future<scalar_data> right_p,
        hpx::shared_future<scalar_data> bottom_p, hpx::shared_future<scalar_data> top_p,
        hpx::shared_future<scalar_data> middle_rhs, 
        std::vector<std::pair<uint, uint> > const& fluid,
        RealType dx_sq, RealType dy_sq, RealType part1, RealType part2);
    
    static RealType compute_residual(
       hpx::shared_future<scalar_data> middle_p,
        hpx::shared_future<scalar_data> left_p, hpx::shared_future<scalar_data> right_p,
        hpx::shared_future<scalar_data> bottom_p, hpx::shared_future<scalar_data> top_p,
        hpx::shared_future<scalar_data> middle_rhs,
        std::vector<std::pair<uint, uint> > const& fluid,
        RealType dx, RealType dy);
    
    static hpx::future<
            std::pair<vector_partition, std::pair<RealType, RealType> >
            >
    update_velocities(
       vector_partition const& middle_uv, scalar_partition const& middle_p,
       scalar_partition const& right_p, scalar_partition const& top_p, 
       vector_partition const& middle_fg,
       std::vector<std::bitset<5> > const& flag_data,
        std::vector<std::pair<uint, uint> > const& fluid,
       RealType dx, RealType dy, RealType dt);
       
    static hpx::future<
            std::tuple<scalar_partition, scalar_partition, scalar_partition>
           >
    compute_stream_vorticity_heat(
        scalar_partition const& middle_stream,
        scalar_partition const& bottom_stream,
        scalar_partition const& middle_vorticity,
        scalar_partition const& middle_heat,
        scalar_partition const& bottom_heat, vector_partition const& middle_uv,
        vector_partition const& right_uv, vector_partition const& top_uv,
        scalar_partition const& middle_temperature,
        scalar_partition const& right_temperature,
        std::vector<std::bitset<5> > const& flag_data, 
        uint global_i, uint global_j, uint i_max, uint j_max, RealType re,
        RealType pr, RealType dx, RealType dy);   
};

}//namespace

#endif
