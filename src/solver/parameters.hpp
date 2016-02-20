#ifndef SOLVER_PARAMETERS_HPP
#define SOLVER_PARAMETERS_HPP

namespace solver {

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
