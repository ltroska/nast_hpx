#ifndef IO_VTK_WRITER_HPP
#define IO_VTK_WRITER_HPP

#include <hpx/include/iostreams.hpp>
#include <hpx/include/runtime.hpp>
#include <hpx/include/util.hpp>

#include "stepper/server/stepper_server.hpp"


namespace io {

typedef std::vector<std::vector<grid::partition_data> > space;

void do_async_write(space u, uint partitions_x, uint partitions_y, uint cells_x, uint cells_y)
{
    std::cout << "p values on locality " << hpx::get_locality_id() << std::endl;
    for (uint j = partitions_y - 1; j <= partitions_y; --j)
    {
        for (uint row = cells_y - 1; row <= cells_y; --row)
        {
            for (uint i = 0; i < partitions_x; ++i)
            {
                for (uint col = 0; col < cells_x; col++)
                {
                    std::cout << u[i][j].get_cell(col, row).p << " ";
                }
            }
            std::cout << std::endl;
        }
    }

    std::cout << std::endl << std::endl;
}

void async_write(space u, uint partitions_x, uint partitions_y, uint cells_x, uint cells_y)
{
    hpx::util::io_service_pool* pool =
        hpx::get_runtime().get_thread_pool("io_pool");

    // ... and schedule the handler to run on one of its OS-threads.
    pool->get_io_service().post(
        hpx::util::bind(&do_async_write, u, partitions_x, partitions_y, cells_x, cells_y));

}

}//namespace io

#endif
