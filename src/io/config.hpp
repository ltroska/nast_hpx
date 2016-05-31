#ifndef NAST_HPX_IO_CONFIG_HPP
#define NAST_HPX_IO_CONFIG_HPP

#include <iostream>
#include <bitset>
#include <vector>

#include "util/typedefs.hpp"
#include "grid/boundary_data.hpp"

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

        Real part1;
        Real part2;

        Real re;
        Real pr;
        Real omega;
        Real tau;
        Real alpha;
        Real beta;
        Real gx;
        Real gy;

        bool vtk;
        Real delta_vec;

        Real t_end;
        Real initial_dt;

        uint sub_iterations;
        uint iter_max;
        Real eps;
        Real eps_sq;

        std::size_t num_localities;
        std::size_t num_localities_x;
        std::size_t num_localities_y;

        std::size_t num_partitions;
        std::size_t num_partitions_x;
        std::size_t num_partitions_y;

        std::size_t num_x_blocks;
        std::size_t num_y_blocks;

        std::size_t cells_x_per_block;
        std::size_t cells_y_per_block;

        std::size_t rank;
        std::size_t idx;
        std::size_t idy;

        grid::boundary_data u_bnd;
        grid::boundary_data v_bnd;
        grid::boundary_data data_type;

        grid::boundary_data temp_bnd;
        grid::boundary_data temp_data_type;
        Real ti;

        bool with_flag_grid;
        std::vector<std::bitset<6> > flag_grid;

        bool with_initial_uv_grid;
        std::vector<std::pair<Real, Real> > initial_uv_grid;

        template <typename Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            ar & i_max & j_max & num_fluid_cells & x_length
                & y_length & dx & dy & over_dx & over_dy
                & dx_sq & dy_sq & part1 & part2 & re & pr & omega & tau & alpha
                & beta & gx & gy & vtk & delta_vec & t_end & initial_dt
                & sub_iterations & iter_max & eps & eps_sq & num_localities
                & num_localities_x & num_localities_y & num_partitions
                & num_partitions_x & num_partitions_y & num_x_blocks
                & num_y_blocks & cells_x_per_block & cells_y_per_block
                & rank & idx & idy & ti & with_flag_grid & with_initial_uv_grid;
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
                << "\n\tnum_localities = " << config.num_localities
                << "\n\tnum_localities_x = " << config.num_localities_x
                << "\n\tnum_localities_y = " << config.num_localities_y
                << "\n\tnum_partitions = " << config.num_partitions
                << "\n\tnum_partitions_x = " << config.num_partitions_x
                << "\n\tnum_partitions_y = " << config.num_partitions_y
                << "\n\tnum_x_blocks = " << config.num_x_blocks
                << "\n\tnum_y_blocks = " << config.num_y_blocks
                << "\n\tnum_cells_x_per_block = " << config.cells_x_per_block
                << "\n\tnum_cells_y_per_block = " << config.cells_y_per_block
                << "\n\trank = " << config.rank
                << "\n\tidx = " << config.idx
                << "\n\tidy = " << config.idy
                << "\n\tnumFluid = " << config.num_fluid_cells
                << "\n\txLength = " << config.x_length
                << "\n\tyLength = " << config.y_length
                << "\nDATA:"
                << "\n\tRe = " << config.re
                << "\n\tPr = " << config.pr
                << "\n\tomega = " << config.omega
                << "\n\tgx = " << config.gx
                << "\n\tgy = " << config.gy
                << "\n\tbeta = " << config.beta
                << "\n\tdt = " << config.initial_dt
                << "\n\tu_bnd " << config.u_bnd
                << "\n\tv_bnd " << config.v_bnd
                << "\n\tbnd_data_type " << config.data_type
                << "\n\ttemp_bnd " << config.temp_bnd
                << "\n\ttemp_data_type " << config.temp_data_type
                << "\n\tt_inital = " << config.ti
                << "\nSIMULATION:"
                << "\n\ttau = " << config.tau
                << "\n\teps = " << config.eps
                << "\n\teps_sq = " << config.eps_sq
                << "\n\talpha = " << config.alpha
                << "\n\tbeta = " << config.beta
                << "\n\tpart1 = " << config.part1
                << "\n\tpart2 = " << config.part2
                << "\n\tdelta_vec = " << config.delta_vec
                << "\n\titer_max = " << config.iter_max
                << "\n\tsub_iterations = " << config.sub_iterations
                << "\n\tvtk = " << config.vtk
                << "\n}";
            return os;
        }

        static config read_config_from_file(const char *path, std::size_t rank, std::size_t num_localities);

};

} //namespace
}

#endif

