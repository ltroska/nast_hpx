#include <hpx/include/iostreams.hpp>

#include "stepper/fd_stepper.hpp"

namespace stepper {

//calculate locality to distribute partitions evenly
inline uint locidx(uint i, uint np, uint nl)
{
    return i / (np/nl);
}

fd_stepper::fd_stepper(cfd_config* config)
{
    this->config = config;
    std::cout << *(this->config);
}

fd_stepper::space fd_stepper::do_work()
{
    std::vector<hpx::id_type> localities = hpx::find_all_localities();
    numLocalities = localities.size();

    uint numPartitionsX = (config->iMax+2)/config->iRes;
    uint numPartitionsY = (config->jMax+2)/config->jRes;
    numPartitions = numPartitionsX * numPartitionsY;

    dx = config->xLength / (config->iMax+1);
    dy = config->yLength / (config->jMax+1);


    std::cout << "#localities = " << numLocalities << " #partitions = " << numPartitions << " dx = " << dx << " dy = " << dy << std::endl;

    space U;

    U.resize(numPartitionsX);

    for (uint i = 0; i < numPartitionsX; ++i)
        U[i].resize(numPartitionsX);

    for (uint j = 0; j < numPartitionsY; ++j)
    {
        for (uint i = 0; i < numPartitionsX; ++i)
        {
            U[i][j] = grid::partition(localities[locidx(j * numPartitionsX + i, numPartitions, numLocalities)],
                                       config->iRes, config->jRes, j*numPartitionsX+i);
        }
    }

    //grid is saved "upside down" (0,0) is in bottom left
    for (uint j = numPartitionsY - 1; j <= numPartitionsY; --j)
    {
        for (uint row = config->jRes-1; row <= config->jRes; --row)
        {
            for (uint i = 0; i < numPartitionsX; ++i)
            {
                for (uint col = 0; col < config->iRes; col++)
                {
                    std::cout << U[i][j].get_data(grid::center_partition).get().get_cell(col, row) << " ";
                }
            }
            std::cout << std::endl;
        }
    }

    std::vector<hpx::future<uint> > set_futures;

    for (uint j = 0; j < numPartitionsY; j++) {
        set_futures.push_back(U[0][j].set_boundary(grid::left_boundary, 4.0));
        set_futures.push_back(U[numPartitionsX-1][j].set_boundary(grid::right_boundary, 2.0));
    }

    for (uint i = 0; i < numPartitionsX; i++) {
        set_futures.push_back(U[i][0].set_boundary(grid::bottom_boundary, 9.0));
        set_futures.push_back(U[i][numPartitionsY-1].set_boundary(grid::top_boundary, 7.0));
    }

    hpx::wait_all(set_futures);

    std::cout << std::endl << std::endl;

    for (uint j = numPartitionsY - 1; j <= numPartitionsY; --j)
    {
        for (uint row = config->jRes-1; row <= config->jRes; --row)
        {
            for (uint i = 0; i < numPartitionsX; ++i)
            {
                for (uint col = 0; col < config->iRes; col++)
                {
                    std::cout << U[i][j].get_data(grid::center_partition).get().get_cell(col, row) << " ";
                }
            }
            std::cout << std::endl;
        }
    }

    std::cout << std::endl << std::endl;

    std::cout << U[3][3].get_data(grid::left_partition).get();
    std::cout << U[3][3].get_data(grid::right_partition).get();
    std::cout << U[3][3].get_data(grid::top_partition).get();
    std::cout << U[3][3].get_data(grid::bottom_partition).get();
}

}//namespace stepper
