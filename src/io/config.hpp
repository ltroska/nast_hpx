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
        RealType pr;
        RealType omega;
        RealType tau;
        RealType eps;
        RealType eps_sq;
        RealType alpha;
        RealType beta;
        RealType gx;
        RealType gy;

        bool vtk;
        uint output_skip_size;
        RealType delta_vec;

        RealType t_end;
        RealType dt;
        uint sub_iterations;
        uint iter_max;

        uint wfe;

        boundary_data u_bnd;
        boundary_data v_bnd;
        boundary_data data_type;

        boundary_data temp_bnd;
        boundary_data temp_data_type;
        RealType ti;

        bool with_flag_grid;
        std::vector<std::bitset<5> > flag_grid;

        bool with_initial_uv_grid;
        std::vector<std::pair<RealType, RealType> > initial_uv_grid;

        template <typename Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            ar & i_max & j_max & i_res & j_res & x_length & y_length & re & pr & omega & tau & eps & eps_sq & alpha & beta
            & t_end & dt & delta_vec & vtk & output_skip_size & sub_iterations & iter_max & wfe & gx & gy; //& u_bnd & v_bnd & data_type
           // & temp_bnd & temp_data_type & ti & with_flag_grid & flag_grid & with_initial_uv_grid & initial_uv_grid;
        }

        friend std::ostream& operator<<(std::ostream& os, config const& config)
        {
            os  << "config:"
                << "\n{"
                << "\n\tiMax = " << config.i_max
                << "\n\tjMax = " << config.j_max
                << "\n\tiRes = " << config.i_res
                << "\n\tjRes = " << config.j_res
                << "\n\txLength = " << config.x_length
                << "\n\tyLength = " << config.y_length
                << "\n\tRe = " << config.re
                << "\n\tPr = " << config.pr
                << "\n\tomega = " << config.omega
                << "\n\ttau = " << config.tau
                << "\n\teps = " << config.eps
                << "\n\talpha = " << config.alpha
                << "\n\toutputSkipSize = " << config.output_skip_size
                << "\n\tdeltaVec = " << config.delta_vec
                << "\n\titerMax = " << config.iter_max
                << "\n\tsubIterations = " << config.sub_iterations
                << "\n\tGX = " << config.gx
                << "\n\tGY = " << config.gy
                << "\n\tbeta = " << config.beta
                << "\n\tdt = " << config.dt
                << "\n\twfe = " << config.wfe
                << "\n\tvtk = " << config.vtk
                << "\n\tu_bnd " << config.u_bnd
                << "\n\tv_bnd " << config.v_bnd
                << " bnd_data_type " << config.data_type
                << "\n\ttemp_bnd " << config.temp_bnd
                << "\n\ttemp_data_type " << config.temp_data_type
                << "\n\tTI = " << config.ti
                << "\n}";
            return os;
        }
};

} //namespace

#endif

