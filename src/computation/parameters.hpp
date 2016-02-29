#ifndef COMPUTATION_PARAMETERS_HPP
#define COMPUTATION_PARAMETERS_HPP

namespace computation {

struct parameters
{
    uint i_max, j_max, num_partitions_x, num_partitions_y, num_cells_per_partition_x, num_cells_per_partition_y;
    RealType re, omega, alpha, dx, dy;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & i_max & j_max & num_partitions_x & num_partitions_y & num_cells_per_partition_x & num_cells_per_partition_y & re & omega & alpha & dx & dy;
    }
};

}//namespace

#endif
