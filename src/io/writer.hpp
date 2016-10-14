#ifndef NAST_HPX_IO_WRITER_HPP_
#define NAST_HPX_IO_WRITER_HPP_

#include "grid/partition_data.hpp"

#include <bitset>

namespace nast_hpx { namespace io {

    typedef grid::partition_data<double> grid_type;
    typedef grid::partition_data<std::bitset<9> > type_grid;

    struct writer
    {
        static void write_vtk(grid_type const& p_data, grid_type const& u_data,
            grid_type const& v_data, grid_type const& w_data, type_grid const& cell_types,
            std::size_t res_x, std::size_t res_y, std::size_t res_z, std::size_t i_max,
            std::size_t j_max, std::size_t k_max, double dx, double dy, double dz, std::size_t step,
            std::size_t loc, std::size_t idx, std::size_t idy, std::size_t idz);
    };

}
}

#endif
