#ifndef NAST_HPX_IO_CONFIG_HPP
#define NAST_HPX_IO_CONFIG_HPP

#include <iostream>
#include <bitset>
#include <vector>

#include "util/defines.hpp"
#include "grid/boundary_condition.hpp"

namespace nast_hpx { namespace io {

/// This class represents the configuration of the simulation provided by the
/// csv file passed in as an argument to the executable.
struct config
{
        uint i_max;
        uint j_max;
        uint num_fluid_cells;

        Real x_length;
        Real y_length;
        Real dx;
        Real dy;
        Real over_dx;
        Real over_dy;
        Real dx_sq;
        Real dy_sq;
        Real over_dx_sq;
        Real over_dy_sq;

        Real re;
        Real omega;
        Real tau;
        Real alpha;
        Real beta;
        Real gx;
        Real gy;

        Real part1;
        Real part2;

        bool vtk;
        Real delta_vec;
        bool verbose;

        Real t_end;
        Real initial_dt;
        std::size_t max_timesteps;

        uint iter_max;
        Real eps;
        Real eps_sq;

        std::size_t num_x_blocks;
        std::size_t num_y_blocks;

        std::size_t cells_x_per_block;
        std::size_t cells_y_per_block;

        grid::boundary_condition bnd_condition;

        std::vector<std::bitset<7> > flag_grid;
        std::vector<std::bitset<5> > empty_marker_grid;

        bool with_initial_uv_grid;
        std::vector<pair<Real> > initial_uv_grid;

        template <typename Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            ar & i_max & j_max & num_fluid_cells & x_length
                & y_length & dx & dy & over_dx & over_dy
                & dx_sq & dy_sq & re & omega & tau & alpha
                & beta & gx & gy & vtk & delta_vec & verbose & t_end & initial_dt & max_timesteps
                & iter_max & eps & eps_sq & part1 & part2
                & num_x_blocks & num_y_blocks & cells_x_per_block & cells_y_per_block
                & with_initial_uv_grid & bnd_condition;
        }

        friend std::ostream& operator<<(std::ostream& os, config const& config)
        {
            os  << "config:"
                << "\n{"
                << "\nGEOMETRY:"
                << "\n\ti_max = " << config.i_max
                << "\n\tj_max = " << config.j_max
                << "\n\tdx = " << config.dx
                << "\n\tdy = " << config.dy
                << "\n\tdx_sq = " << config.dx_sq
                << "\n\tdy_sq = " << config.dy_sq
                << "\n\tover_dx = " << config.over_dx
                << "\n\tover_dy = " << config.over_dy
                << "\n\tover_dx_sq = " << config.over_dx_sq
                << "\n\tover_dy_sq = " << config.over_dy_sq
                << "\n\tnum_x_blocks = " << config.num_x_blocks
                << "\n\tnum_y_blocks = " << config.num_y_blocks
                << "\n\tnum_cells_x_per_block = " << config.cells_x_per_block
                << "\n\tnum_cells_y_per_block = " << config.cells_y_per_block
                << "\n\tnumFluid = " << config.num_fluid_cells
                << "\n\txLength = " << config.x_length
                << "\n\tyLength = " << config.y_length
                << "\nDATA:"
                << "\n\tRe = " << config.re
                << "\n\tomega = " << config.omega
                << "\n\tgx = " << config.gx
                << "\n\tgy = " << config.gy
                << "\n\tbeta = " << config.beta
                << "\n\tpart1 = " << config.part1
                << "\n\tpart2 = " << config.part2
                << "\n\tdt = " << config.initial_dt
                << "\n\tt_end = " << config.t_end
                << "\n\tmax_timesteps = " << config.max_timesteps
                << "\n\tboundary_condition = " << config.bnd_condition
                << "\nSIMULATION:"
                << "\n\ttau = " << config.tau
                << "\n\teps = " << config.eps
                << "\n\teps_sq = " << config.eps_sq
                << "\n\talpha = " << config.alpha
                << "\n\tbeta = " << config.beta
                << "\n\tdelta_vec = " << config.delta_vec
                << "\n\titer_max = " << config.iter_max
                << "\n\tvtk = " << config.vtk
                << "\n}";
            return os;
        }

        static config read_config_from_file(const char *path, std::size_t rank, std::size_t num_localities);

};

} //namespace
}

#endif

