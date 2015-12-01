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

    dx = config->xLength / (config->iMax + 1);
    dy = config->yLength / (config->jMax + 1);


    std::cout << "#localities = " << numLocalities << " #partitions = " << numPartitions << " dx = " << dx << " dy = " << dy << std::endl;

    space U;

    grid::partition a = grid::partition(hpx::find_here(), config->iRes, config->jRes, 1);

    std::cout << a.get_data(grid::center_partition).get() << std::endl;

  /*  for (uint j = 0; j < numPartitionsY; ++j)
    {
        for (uint i = 0; i < numPartitionsX; ++i)
        {
            std::cout << localities[locidx(j * numPartitionsX + i, numPartitions, numLocalities)] << std::endl;
            U[i][j] = grid::partition(hpx::find_here(),
                                       config->iRes, config->jRes, 1);
        }
    }

    std::cout << "HAHA" << std::endl << hpx::flush;


    for (uint j = numPartitionsY - 1; j != 0; --j)
    {
        for (uint i = 0; i < numPartitionsX; ++i)
        {
            std::cout << "HAHA" << std::endl;
            std::cout << U[i][j].get_data(grid::center_partition).get() << std::endl;
        }
    }*/
}

}//namespace stepper
