#ifndef TEST_HELPERS_HPP
#define TEST_HELPERS_HPP

#include "grid/partition.hpp"

/*
* Makes the needed grids for testing the computation component
*/
struct grid_maker
{
    uint num_localities_x, num_localities_y;
    uint locality_id;

    computation::parameters p;

    grid_maker(uint num_localities_x_, uint num_localities_y_, uint locality_id_, computation::parameters params)
            : num_localities_x(num_localities_x_), num_localities_y(num_localities_y_), locality_id(locality_id_), p(params)
    {}

    void make_index_grid(index_grid_type& grid)
    {
        grid.resize(p.num_partitions_x * p.num_partitions_y);

        for (uint l = 0; l < p.num_partitions_y; l++)
            for (uint k = 0; k < p.num_partitions_x; k++)
            {
                grid[get_index(k, l)] =
                    std::pair<RealType, RealType>
                        (
                        (locality_id % num_localities_x) * (p.num_partitions_x - 2) * p.num_cells_per_partition_x + (k - 1) * p.num_cells_per_partition_x,
                        (locality_id / num_localities_x) * (p.num_partitions_y - 2) * p.num_cells_per_partition_y + (l - 1) * p.num_cells_per_partition_y
                        );
            }
    }

    void make_vector_grid(vector_grid_type& grid, RealType initial_value)
    {
        grid.resize(p.num_partitions_x * p.num_partitions_y);

        for (uint l = 0; l < p.num_partitions_y; l++)
            for (uint k = 0; k < p.num_partitions_x; k++)
            {
                grid[get_index(k, l)] = vector_partition(hpx::find_here(), p.num_cells_per_partition_x, p.num_cells_per_partition_y, initial_value);
            }
    }

    private:
        uint get_index(uint k, uint l) { return l * p.num_partitions_x + k;}
};

#endif
