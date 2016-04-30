#ifndef CELL_OPERATIONS_HPP
#define CELL_OPERATIONS_HPP

#include <hpx/hpx.hpp>

#include "grid/types.hpp"
#include "util/typedefs.hpp"
#include "util/boundary_data.hpp"


namespace computation
{    
void set_velocity_for_cell(vector_cell& middle,
        vector_cell const& left, vector_cell const& right,
        vector_cell const& bottom, vector_cell const& top,
        std::bitset<5> const& cell_type,
        boundary_data const& type,
        boundary_data const& u,
        boundary_data const& v);

void set_temperature_for_cell(scalar_cell& middle,
        scalar_cell const& left, scalar_cell const& right,
        scalar_cell const& bottom, scalar_cell const& top,
        boundary_data const& type,
        boundary_data const& temperature,
        std::bitset<5> const& cell_type,
        uint i, uint j,
        RealType dx, RealType dy);

void compute_fg_for_cell(vector_cell& middle_fg,
    vector_cell const& middle_uv, vector_cell const& left_uv,
    vector_cell const& right_uv, vector_cell const& bottom_uv,
    vector_cell const& top_uv, vector_cell const& bottomright_uv,
    vector_cell const& topleft_uv, scalar_cell const& middle_temperature,
    scalar_cell const& right_temperature, scalar_cell const& top_temperature,
    std::bitset<5> const& type,
    RealType re, RealType gx, RealType gy, RealType beta,
    RealType dx, RealType dy, RealType dt, RealType alpha);

void compute_temperature_for_cell(scalar_cell& middle_temperature,
    scalar_cell const& old_middle_temperature,
    scalar_cell const& left_temperature, scalar_cell const& right_temperature,
    scalar_cell const& bottom_temperature, scalar_cell const& top_temperature,
    vector_cell const& middle_uv, vector_cell const& left_uv,
    vector_cell const& bottom_uv, std::bitset<5> const& type,
    RealType re, RealType pr, RealType dx, RealType dy, RealType dt,
    RealType alpha);

void compute_rhs_for_cell(scalar_cell& middle_rhs,
    vector_cell const& middle_fg, vector_cell const& left_fg,
    vector_cell const& bottom_fg, std::bitset<5> const& type,
    RealType dx, RealType dy, RealType dt);

void set_pressure_for_cell(scalar_cell& middle_p,
    scalar_cell const& left_p, scalar_cell const& right_p,
    scalar_cell const& bottom_p, scalar_cell const& top_p,
    std::bitset<5> const& type);

void do_sor_cycle_for_cell(scalar_cell& middle_p,
    scalar_cell const& left_p, scalar_cell const& right_p,
    scalar_cell const& bottom_p, scalar_cell const& top_p,
    scalar_cell const& middle_rhs, std::bitset<5> const& type,
    RealType dx_sq, RealType dy_sq, RealType part1, RealType part2);

RealType compute_residual_for_cell(scalar_cell const& middle_p,
    scalar_cell const& left_p, scalar_cell const& right_p,
    scalar_cell const& bottom_p, scalar_cell const& top_p,
    scalar_cell const& middle_rhs, std::bitset<5> const& type,
    RealType over_dx_sq, RealType over_dy_sq);

void update_velocity_for_cell(vector_cell& middle_uv,
    scalar_cell const& middle_p, scalar_cell const& right_p,
    scalar_cell const& top_p, vector_cell const& middle_fg,
    std::bitset<5> const& type, RealType over_dx, RealType over_dy,
    RealType dt);
}//end namespace computation

#endif /* CELL_OPERATIONS_HPP */

