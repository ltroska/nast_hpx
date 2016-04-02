#ifndef COMPUTATION_STRATEGY_HPP
#define COMPUTATION_STRATEGY_HPP

#include <hpx/hpx.hpp>

#include "grid/types.hpp"
#include "util/boundary_data.hpp"

namespace computation {

class strategy {
    public:
        virtual void set_boundary(vector_data& uv_center, vector_data const& uv_left, vector_data const& uv_right, vector_data const& uv_bottom, vector_data const& uv_top,
                                    scalar_data& temperature, std::vector<std::bitset<5> > const& flag_data, boundary_data const& data_type, boundary_data const& u_bnd,
                                    boundary_data const& v_bnd, boundary_data const& temp_bnd, uint global_i, uint global_j, uint i_max, uint j_max) {}

        virtual void compute_fg(vector_data& fg, vector_data const& uv_center,
                                    vector_data const& uv_left, vector_data const& uv_right,
                                    vector_data const& uv_bottom, vector_data const& uv_top,
                                    vector_data const& uv_bottomright, vector_data const& uv_topleft,
                                    std::vector<std::bitset<5> > const& flag_data,
                                    uint global_i, uint global_j, uint i_max, uint j_max, RealType re,
                                    RealType dx, RealType dy, RealType dt, RealType alpha) {}

        virtual void compute_rhs(scalar_data& rhs, vector_data const& fg_center, vector_data const& fg_left,
                                    vector_data const& fg_bottom, std::vector<std::bitset<5> > const& flag_data, uint global_i, uint global_j, uint i_max, uint j_max,
                                    RealType dx, RealType dy, RealType dt) {}

        virtual void set_pressure_on_boundary(scalar_data& p_center, scalar_data const& p_left, scalar_data const& p_right,
                                                scalar_data const& p_bottom, scalar_data const& p_top, std::vector<std::bitset<5> > const& flag_data,
                                                uint global_i, uint global_j, uint i_max, uint j_max) {}

        virtual void sor_cycle(scalar_data& p_center, scalar_data const& p_left, scalar_data const& p_right,
                                    scalar_data const& p_bottom, scalar_data const& p_top,
                                    scalar_data const& rhs_center, std::vector<std::bitset<5> > const& flag_data,
                                    uint global_i, uint global_j, uint i_max, uint j_max,
                                    RealType omega, RealType dx, RealType dy) {}

        virtual RealType compute_residual(scalar_data const& p_center, scalar_data const& p_left,
                                            scalar_data const& p_right, scalar_data const& p_bottom,
                                            scalar_data const& p_top, scalar_data const& rhs_center,
                                            std::vector<std::bitset<5> > const& flag_data,
                                            uint global_i, uint global_j, uint i_max, uint j_max, RealType dx,
                                            RealType dy) {}

        virtual void update_velocities(vector_data& uv_center, scalar_data const& p_center, scalar_data const& p_right,
                                        scalar_data const& p_top, vector_data const& fg_center, std::vector<std::bitset<5> > const& flag_data,
                                        uint global_i, uint global_j, uint i_max, uint j_max, RealType dx, RealType dy, RealType dt) {}

        virtual void compute_stream_and_vorticity(scalar_data& stream_center, scalar_data& vorticity_center, scalar_data const& stream_bottom,
                                                    vector_data const& uv_center, vector_data const& uv_right, vector_data const& uv_top,
                                                    uint global_i, uint global_j, uint i_max, uint j_max, RealType dx, RealType dy) {}



};

}//namespace

#endif
