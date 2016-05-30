#ifndef NAST_HPX_IO_WRITER_HPP_
#define NAST_HPX_IO_WRITER_HPP_

#include "grid/partition_data.hpp"

namespace nast_hpx { namespace io {
 
    typedef grid::partition_data<Real> grid_type;
    
    struct writer
    {
        static void write_vtk(grid_type const& p_data, grid_type const& u_data, grid_type const& v_data,
            std::size_t res_x, std::size_t res_y, std::size_t i_max,
            std::size_t j_max, Real dx, Real dy, std::size_t step, std::size_t loc);
    };
   
}
}

#endif