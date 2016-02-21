#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include "grid/partition.hpp"
#include "computation/custom_grain_size.hpp"
#include "util/cell.hpp"
#include "test_helpers.hpp"

bool check_set_boundary(vector_grid_type grid, index_grid_type index, computation::parameters p)
{
    std::vector<std::vector<vector_data > > data;

    data.resize(p.num_partitions_x);
    for (uint k = 0; k < p.num_partitions_x ; k++)
    {
        data[k].resize(p.num_partitions_y);
        for (uint l = 0; l < p.num_partitions_y ; l++)
        {
            vector_data base = grid[l * p.num_partitions_x + k].get_data(CENTER).get();
            data[k][l] = vector_data(base);
        }
    }

    uint global_i, global_j, size_x, size_y;

    for (uint l = 1; l < p.num_partitions_y-1; l++)
    {
        for (uint k = 1; k < p.num_partitions_x-1; k++)
        {
            size_x = data[k][l].size_x();
            size_y = data[k][l].size_y();

            for (uint j = 0; j < size_y; j++)
            {
                for (uint i = 0; i < size_x; i++)
                {
                    global_i = index[l * p.num_partitions_x + k].first + i;
                    global_j = index[l * p.num_partitions_x + k].second + j;

                    vector_cell cell = data[k][l].get_cell(i, j);

                    if (global_i == 0 && (global_j >= 1 && global_j <= p.j_max))
                    {
                        if (cell.first != 0)
                            return false;
                        if (cell.second != -data[k][l].get_cell(i+1, j).second)
                            return false;
                    }

                    if (global_i == p.i_max && (global_j >= 1 && global_j <= p.j_max))
                    {
                        if (cell.first != 0)
                            return false;
                        if (cell.second != -data[k][l].get_cell(i+1, j).second)
                            return false;
                    }

                    if (global_j == 0 && (global_i >= 1 && global_i <= p.i_max))
                    {
                        if (cell.second != 0)
                            return false;
                        if (cell.first != -data[k][l].get_cell(i, j+1).first)
                            return false;
                    }

                    if (global_j == p.j_max && (global_i >= 1 && global_i <= p.i_max))
                    {
                        if (cell.second != 0)
                            return false;
                        if (cell.first != 2-data[k][l].get_cell(i, j+1).first)
                            return false;
                    }

                }
            }
        }
    }

    return true;

}

void do_uv_test(uint i_max, uint j_max, uint locality_id, uint localities_x, uint localities_y, uint i_res, uint j_res)
{
        uint num_partitions_x = i_res + 2;
        uint num_partitions_y = j_res + 2;

        computation::parameters params;
        params.num_partitions_x = num_partitions_x;
        params.num_partitions_y = num_partitions_y;
        params.i_max = i_max;
        params.j_max = j_max;
        params.num_cells_per_partition_x = (i_max / localities_x) / i_res;
        params.num_cells_per_partition_y = (j_max / localities_y) / j_res;

        grid_maker maker = grid_maker(localities_x, localities_y, locality_id, params);

        index_grid_type index;
        vector_grid_type grid;

        maker.make_index_grid(index);
        maker.make_vector_grid(grid, 0);

        computation::custom_grain_size strat(index, params);
        strat.set_velocity_on_boundary(grid);

        std::string msg = "\nfailed with settings " + std::to_string(i_max) + " " + std::to_string(j_max) + " " + std::to_string(locality_id) + " "
                            + std::to_string(localities_x) + " " + std::to_string(localities_y) + " " + std::to_string(i_res) + " " + std::to_string(j_res) + "\n";

        HPX_ASSERT_MSG(check_set_boundary(grid, index, params), msg.c_str());
}

int hpx_main(int argc, char* argv[])
{
// ------------------------------- SET VELOCITY ON BOUNDARY ------------------------------- //

// --- SQUARE AREA, ONE PARTITION --- //
    do_uv_test(20, 20, 0, 1, 1, 1, 1);
    do_uv_test(64, 64, 0, 1, 1, 1, 1);

// --- SQUARE AREA, MULTIPLE PARTITIONS --- //
    do_uv_test(64, 64, 0, 1, 1, 8, 8);
    do_uv_test(256, 256, 0, 1, 1, 8, 16);
    do_uv_test(256, 256, 0, 1, 1, 1, 8);
    do_uv_test(256, 256, 0, 1, 1, 4, 8);

// --- SQUARE AREA, MULTIPLE PARTITIONS, NOT FIRST LOCALITY --- //
    do_uv_test(32, 32, 2, 2, 2, 8, 8);
    do_uv_test(1028, 1028, 1, 1, 2, 8, 16);
    do_uv_test(1028, 1028, 1, 2, 1, 8, 16);
    do_uv_test(2056, 2056, 3, 2, 2, 1, 8);
    do_uv_test(2056, 2056, 7, 4, 4, 4, 8);

// --- NON-SQUARE AREA, ONE PARTITION --- //
    do_uv_test(20, 40, 0, 1, 1, 1, 1);
    do_uv_test(40, 20, 0, 1, 1, 1, 1);
    do_uv_test(64, 512, 0, 1, 1, 1, 1);
    do_uv_test(512, 64, 0, 1, 1, 1, 1);

// --- NON-SQUARE AREA, MULTIPLE PARTITIONS --- //
    do_uv_test(64, 32, 0, 1, 1, 8, 8);
    do_uv_test(32, 64, 0, 1, 1, 8, 8);
    do_uv_test(128, 256, 0, 1, 1, 8, 16);
    do_uv_test(256, 128, 0, 1, 1, 8, 16);
    do_uv_test(256, 128, 0, 1, 1, 1, 16);
    do_uv_test(256, 128, 0, 1, 1, 1, 16);
    do_uv_test(256, 64, 0, 1, 1, 4, 8);
    do_uv_test(64, 256, 0, 1, 1, 4, 8);

// --- NON-SQUARE AREA, MULTIPLE PARTITIONS, NOT FIRST LOCALITY --- //
    do_uv_test(64, 32, 2, 2, 2, 8, 8);
    do_uv_test(32, 64, 2, 2, 2, 8, 8);
    do_uv_test(1028, 512, 1, 1, 2, 8, 16);
    do_uv_test(512, 1028, 1, 1, 2, 8, 16);
    do_uv_test(1028, 2056, 1, 2, 1, 8, 16);
    do_uv_test(2056, 1028, 1, 2, 1, 8, 16);
    do_uv_test(2056, 512, 3, 2, 2, 1, 8);
    do_uv_test(512, 2056, 3, 2, 2, 1, 8);
    do_uv_test(2056, 1028, 7, 4, 4, 4, 8);
    do_uv_test(1028, 2056, 7, 4, 4, 4, 8);

// ---------------------------------------------------------------------------------------- //

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}
