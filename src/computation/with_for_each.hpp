/** Strategy implements the algorithm using the parallel algorithms already
 * implemented in HPX.
 */

#ifndef WITH_FOR_EACH_HPP
#define WITH_FOR_EACH_HPP

#include <hpx/hpx.hpp>

#include "grid/types.hpp"
#include "util/boundary_data.hpp"

namespace computation {

class with_for_each {

public:
    static vector_partition set_velocity_for_boundary_and_obstacles(
            vector_partition const& middle_partition,
            vector_partition const& left_partition,
            vector_partition const& right_partition,
            vector_partition const& bottom_partition,
            vector_partition const& top_partition,
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
            std::vector<std::bitset<5> > const& flag_data,
            boundary_data const& boundary_data_type,
            boundary_data const& temperature_boundary_data,
            uint global_i, uint global_j,
            Real dx, Real dy);
    
    static vector_partition compute_fg_on_fluid_cells(
        vector_partition const& middle_uv,
        vector_partition const& left_uv, vector_partition const& right_uv,
        vector_partition const& bottom_uv, vector_partition const& top_uv,
        vector_partition const& bottomright_uv,
        vector_partition const& topleft_uv,
        scalar_partition const& middle_temperature,
        scalar_partition const& right_temperature,
        scalar_partition const& top_temperature,
        std::vector<std::bitset<5> > const& flag_data,
        Real re, Real gx, Real gy, Real beta,
        Real dx, Real dy, Real dt, Real alpha);
    
    static scalar_partition compute_temperature_on_fluid_cells(
        scalar_partition const& middle_temperature,
        scalar_partition const& left_temperature,
        scalar_partition const& right_temperature,
        scalar_partition const& bottom_temperature,
        scalar_partition const& top_temperature,
        vector_partition const& middle_uv,
        vector_partition const& left_uv,
        vector_partition const& bottom_uv,
        std::vector<std::bitset<5> > const& flag_data,
        Real re, Real pr, Real dx, Real dy, Real dt,
        Real alpha);
    
    static scalar_partition compute_right_hand_side_on_fluid_cells(
        vector_partition const& middle_fg, vector_partition const& left_fg,
        vector_partition const& bottom_fg,
        std::vector<std::bitset<5> > const& flag_data,
        Real dx, Real dy, Real dt);
    
    static scalar_partition set_pressure_for_boundary_and_obstacles(
        scalar_partition const& middle_p, scalar_partition const& left_p,
        scalar_partition const& right_p, scalar_partition const& bottom_p,
        scalar_partition const& top_p,
        std::vector<std::bitset<5> > const& flag_data);
    
    static scalar_partition sor_cycle(scalar_partition const& middle_p,
        scalar_partition const& left_p, scalar_partition const& right_p,
        scalar_partition const& bottom_p, scalar_partition const& top_p,
        scalar_partition const& middle_rhs, 
        std::vector<std::bitset<5> > const& flag_data,
        Real omega, Real dx, Real dy);
    
  
    static hpx::future<Real> compute_residual(
        scalar_partition const& middle_p, scalar_partition const& left_p,
        scalar_partition const& right_p, scalar_partition const& bottom_p,
        scalar_partition const& top_p, scalar_partition const& middle_rhs,
        std::vector<std::bitset<5> > const& flag_data, Real dx,
        Real dy);
    
    static hpx::future<
            std::pair<vector_partition, std::pair<Real, Real> >
            >
    update_velocities(
       vector_partition const& middle_uv, scalar_partition const& middle_p,
       scalar_partition const& right_p, scalar_partition const& top_p, 
       vector_partition const& middle_fg,
       std::vector<std::bitset<5> > const& flag_data,
       Real dx, Real dy, Real dt);
       
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
        uint global_i, uint global_j, uint i_max, uint j_max, Real re,
        Real pr, Real dx, Real dy);   
};

}//namespace

#endif
