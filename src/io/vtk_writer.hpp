#ifndef IO_VTK_WRITER_HPP
#define IO_VTK_WRITER_HPP

#include <hpx/include/iostreams.hpp>
#include <hpx/include/runtime.hpp>
#include <hpx/include/util.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/thread_executors.hpp>

#include <fstream>

#include "stepper/server/stepper_server.hpp"

namespace io {

template<typename T>
void do_async_print(std::vector<std::vector<grid::partition_data<T> > > const& grid, std::string const message, uint partitions_x, uint partitions_y, uint cells_x, uint cells_y, boost::shared_ptr<hpx::lcos::local::promise<int> > p)
{
    uint res_x_, res_y_;
    if (hpx::get_num_localities_sync() == 2)
    {
        res_x_ = 2;
        res_y_ = 1;
    }
    else
    {
        res_x_ = static_cast<uint>(sqrt(hpx::get_num_localities_sync()));
        res_y_ = res_x_;
    }

    std::cout << message  << std::endl;
    for (uint j = partitions_y - 1; j <= partitions_y; --j)
    {
        for (uint row = cells_y - 1 ; row <= cells_y; --row)
        {
            for (uint i = 0; i < partitions_x; ++i)
            {
                for (uint col = 0; col <  cells_x; col++)
                {
                    std::cout << grid[i][j].get_cell(col, row) << " ";
                }
            }
            std::cout << std::endl;
        }
    }

    p->set_value(0);
}
}//namespace io

#endif
