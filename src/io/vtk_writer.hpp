#ifndef IO_VTK_WRITER_HPP
#define IO_VTK_WRITER_HPP

#include <hpx/include/iostreams.hpp>
#include <hpx/include/runtime.hpp>
#include <hpx/include/util.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/thread_executors.hpp>

#include "stepper/server/stepper_server.hpp"

namespace io {

typedef std::vector<std::vector<grid::partition_data> > space;

void do_async_write(space u, uint partitions_x, uint partitions_y, uint cells_x, uint cells_y, boost::shared_ptr<hpx::lcos::local::promise<int> > p)
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

    std::cout << "p values on locality " << hpx::get_locality_id() << std::endl;
    for (uint j = partitions_y - 1; j <= partitions_y; --j)
    //for (uint j = partitions_y - 2; j > 0; --j)
    {
        for (uint row = ((j == 0 || j == partitions_y -1)) ? 0 : cells_y - 1 ; row <= ((j == 0 || j == partitions_y -1) ? 0 : cells_y - 1); --row)
        {
            for (uint i = 0; i < partitions_x; ++i)
            //for (uint i = 1; i < partitions_x - 1; ++i)
            {
                for (uint col = 0; col < ((i == 0 || i == partitions_x -1 ) ? 1 : cells_x); col++)
                {
                    std::cout << u[i][j].get_cell(col, row).p << " ";
                }
            }
            std::cout << std::endl;
        }
    }

    std::cout << std::endl << std::endl;

    std::cout << "uv values on locality " << hpx::get_locality_id() << std::endl;
    for (uint j = partitions_y - 1; j <= partitions_y; --j)
    //for (uint j = partitions_y - 2; j > 0; --j)
    {
        for (uint row = ((j == 0 || j == partitions_y -1)) ? 0 : cells_y - 1 ; row <= ((j == 0 || j == partitions_y -1) ? 0 : cells_y - 1); --row)
        {
            for (uint i = 0; i < partitions_x; ++i)
            //for (uint i = 1; i < partitions_x - 1; ++i)
            {
                for (uint col = 0; col < ((i == 0 || i == partitions_x -1 ) ? 1 : cells_x); col++)
                {
                    std::cout << "{" << u[i][j].get_cell(col, row).u << "|" << u[i][j].get_cell(col, row).v << "} ";
                }
            }
            std::cout << std::endl;
        }
    }

    std::cout << std::endl << std::endl;

        std::cout << std::endl << std::endl;

    std::cout << "fg values on locality " << hpx::get_locality_id() << std::endl;
    for (uint j = partitions_y - 1; j <= partitions_y; --j)
    //for (uint j = partitions_y - 2; j > 0; --j)
    {
        for (uint row = ((j == 0 || j == partitions_y -1)) ? 0 : cells_y - 1 ; row <= ((j == 0 || j == partitions_y -1) ? 0 : cells_y - 1); --row)
        {
            for (uint i = 0; i < partitions_x; ++i)
            //for (uint i = 1; i < partitions_x - 1; ++i)
            {
                for (uint col = 0; col < ((i == 0 || i == partitions_x -1 ) ? 1 : cells_x); col++)
                {
                    std::cout << "{" << u[i][j].get_cell(col, row).f << "|" << u[i][j].get_cell(col, row).g << "} ";
                }
            }
            std::cout << std::endl;
        }
    }

    std::cout << std::endl << std::endl;

    p->set_value(0);
}

hpx::future<int> async_write(space u, uint partitions_x, uint partitions_y, uint cells_x, uint cells_y)
{
    boost::shared_ptr<hpx::lcos::local::promise<int> > p =
            boost::make_shared<hpx::lcos::local::promise<int> >();

    hpx::util::io_service_pool* pool =
        hpx::get_runtime().get_thread_pool("io_pool");

        // ... and schedule the handler to run on one of its OS-threads.
    pool->get_io_service().post(hpx::util::bind(&do_async_write, u, partitions_x, partitions_y, cells_x, cells_y, p));

    return p->get_future();


}

}//namespace io

#endif
