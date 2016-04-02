#pragma once
#ifndef IO_CONFIG_HPP
#define IO_CONFIG_HPP

#include <bitset>

#include "util/typedefs.hpp"
#include "util/boundary_data.hpp"

namespace io {

struct config
{
        uint i_max;
        uint j_max;
        uint i_res;
        uint j_res;
        RealType x_length;
        RealType y_length;
        RealType re;
        RealType omega;
        RealType tau;
        RealType eps;
        RealType eps_sq;
        RealType alpha;
        RealType t_end;
        RealType dt;
        RealType delta_vec;
        bool vtk;
        uint output_skip_size;
        uint sub_iterations;
        uint iter_max;

        uint wfe;

        RealType gx;
        RealType gy;

        boundary_data u_bnd;
        boundary_data v_bnd;
        boundary_data data_type;
        boundary_data temp_bnd;

        bool with_flag_grid;
        std::vector<std::bitset<5> > flag_grid;

        template <typename Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            ar & i_max & j_max & i_res & j_res & x_length & y_length & re & omega & tau & eps & alpha & t_end & dt & output_skip_size & iter_max;
        }

        friend std::ostream& operator<<(std::ostream& os, config const& config)
        {
            os << "config: {iMax = " << config.i_max << ", jMax = "<< config.j_max
                << ", iRes = " << config.i_res << ", jRes = " << config.j_res << ", xLength = " << config.x_length
                << ", yLength = " << config.y_length << ", Re = " << config.re << ", omega = " << config.omega
                << ", tau = " << config.tau << " eps = " << config.eps << ", alpha = " << config.alpha
                << ", outputSkipSize = " << config.output_skip_size << ", deltaVec = " << config.delta_vec << ", iterMax = " << config.iter_max
                << ", subIterations = " << config.sub_iterations << ", GX = " << config.gx << ", GY = " << config.gy
                << ", dt = " << config.dt << ", wfe = " << config.wfe << ", vtk = " << config.vtk << ", u_bnd " << config.u_bnd << ", v_bnd " << config.v_bnd
                << ", temp_bnd " << config.temp_bnd << " bnd_data_type " << config.data_type << "}";
            return os;
        }
};

} //namespace

#endif

