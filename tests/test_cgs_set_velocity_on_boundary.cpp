#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include "grid/partition.hpp"
#include "computation/custom_grain_size.hpp"
#include "util/cell.hpp"
#include "test_helpers.hpp"

void check_set_velocity_on_boundary(vector_grid_type const& orig_grid, vector_grid_type const& grid, index_grid_type const& index,
                                        computation::parameters p, std::string msg)
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

    std::vector<std::vector<vector_data > > orig_data;

    orig_data.resize(p.num_partitions_x);
    for (uint k = 0; k < p.num_partitions_x ; k++)
    {
        orig_data[k].resize(p.num_partitions_y);
        for (uint l = 0; l < p.num_partitions_y ; l++)
        {
            vector_data base = orig_grid[l * p.num_partitions_x + k].get_data(CENTER).get();
            orig_data[k][l] = vector_data(base);
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
                    std::string ident = " partition/cell " + std::to_string(k) + " " + std::to_string(l) + " " + std::to_string(i) + " " + std::to_string(j) + " ";

                    if (in_range(0, 0, 1, p.j_max, global_i, global_j))
                    {
                        HPX_ASSERT_MSG(cell.first == 0, (msg + ident + expected_string(0, cell.first)).c_str());
                        HPX_ASSERT_MSG(cell.second == -orig_data[k][l].get_cell(i+1, j).second,
                                        (msg + ident + expected_string(-orig_data[k][l].get_cell(i+1, j).second, cell.second)).c_str());
                    }

                    if (in_range(p.i_max, p.i_max, 1, p.j_max, global_i, global_j))
                    {
                        HPX_ASSERT_MSG(cell.first == 0, (msg + ident + expected_string(0, cell.first)).c_str());
                    }

                    if (in_range(p.i_max+1, p.i_max+1, 1, p.j_max, global_i, global_j))
                    {
                        HPX_ASSERT_MSG(cell.second == -orig_data[k][l].get_cell(i-1, j).second,
                                        (msg + ident + expected_string(-orig_data[k][l].get_cell(i-1, j).second, cell.second)).c_str());
                    }

                    if (in_range(1, p.i_max, 0, 0, global_i, global_j))
                    {
                        HPX_ASSERT_MSG(cell.second == 0, (msg + ident + expected_string(0, cell.second)).c_str());
                        HPX_ASSERT_MSG(cell.first == -orig_data[k][l].get_cell(i, j+1).first,
                                        (msg + ident + expected_string(-orig_data[k][l].get_cell(i, j+1).first, cell.first)).c_str());
                    }

                    if (in_range(1, p.i_max, p.j_max, p.j_max, global_i, global_j))
                    {
                        HPX_ASSERT_MSG(cell.second == 0, (msg + ident + expected_string(0, cell.second)).c_str());
                    }

                    if (in_range(1, p.i_max, p.j_max+1, p.j_max+1, global_i, global_j))
                    {
                        HPX_ASSERT_MSG(cell.first == 2-orig_data[k][l].get_cell(i, j-1).first,
                                        (msg + ident + expected_string(2-orig_data[k][l].get_cell(i, j-1).first,cell.first)).c_str());
                    }

                }
            }
        }
    }
}

void do_set_velocity_on_boundary_test(uint i_max, uint j_max, uint locality_id, uint localities_x, uint localities_y, uint i_res, uint j_res)
{

        computation::parameters params;
        params.num_cells_per_partition_x = i_res;
        params.num_cells_per_partition_y = j_res;

        params.num_partitions_x = ((i_max + 2) / localities_x) / i_res + 2;
        params.num_partitions_y = ((j_max + 2) / localities_y) / j_res + 2;

        params.i_max = i_max;
        params.j_max = j_max;

        grid_maker maker = grid_maker(localities_x, localities_y, locality_id, params);

        index_grid_type index;
        vector_grid_type grid, orig_grid;

        maker.make_index_grid(index);
        maker.make_random_grid(grid);

       // print_grid(grid, params);

        maker.copy_grid(grid, orig_grid);

        computation::custom_grain_size strat(index, params);
        strat.set_velocity_on_boundary(grid);

        std::string msg = "\nfailed with settings " + std::to_string(i_max) + " " + std::to_string(j_max) + " " + std::to_string(locality_id) + " "
                            + std::to_string(localities_x) + " " + std::to_string(localities_y) + " " + std::to_string(i_res) + " " + std::to_string(j_res) + " ";

       check_set_velocity_on_boundary(orig_grid, grid, index, params, msg);
}

int hpx_main(int argc, char* argv[])
{
// --- SQUARE AREA, ONE PARTITION --- //
    do_set_velocity_on_boundary_test(6, 6, 0, 1, 1, 8, 8);
    do_set_velocity_on_boundary_test(62, 62, 0, 1, 1, 64, 64);

// --- SQUARE AREA, MULTIPLE PARTITIONS --- //
    do_set_velocity_on_boundary_test(62, 62, 0, 1, 1, 8, 8);
    do_set_velocity_on_boundary_test(254, 254, 0, 1, 1, 8, 16);
    do_set_velocity_on_boundary_test(254, 254, 0, 1, 1, 2, 8);
    do_set_velocity_on_boundary_test(254, 254, 0, 1, 1, 4, 8);

// --- SQUARE AREA, MULTIPLE PARTITIONS, NOT FIRST LOCALITY --- //
    do_set_velocity_on_boundary_test(30, 30, 2, 2, 2, 8, 8);
    do_set_velocity_on_boundary_test(510, 510, 1, 1, 2, 8, 16);
    do_set_velocity_on_boundary_test(510, 510, 1, 2, 1, 8, 16);
   // do_set_velocity_on_boundary_test(1022, 1022, 3, 2, 2, 2, 8);
   // do_set_velocity_on_boundary_test(1022, 1022, 7, 4, 4, 4, 8);

/*
// --- NON-SQUARE AREA, ONE PARTITION --- //
    do_set_velocity_on_boundary_test(14, 30, 0, 1, 1, 1, 1);
    do_set_velocity_on_boundary_test(30, 14, 0, 1, 1, 1, 1);
    do_set_velocity_on_boundary_test(62, 510, 0, 1, 1, 1, 1);
    do_set_velocity_on_boundary_test(510, 62, 0, 1, 1, 1, 1);

// --- NON-SQUARE AREA, MULTIPLE PARTITIONS --- //
    do_set_velocity_on_boundary_test(62, 30, 0, 1, 1, 8, 8);
    do_set_velocity_on_boundary_test(30, 62, 0, 1, 1, 8, 8);
    do_set_velocity_on_boundary_test(126, 254, 0, 1, 1, 8, 16);
    do_set_velocity_on_boundary_test(254, 126, 0, 1, 1, 8, 16);
    do_set_velocity_on_boundary_test(254, 126, 0, 1, 1, 1, 16);
    do_set_velocity_on_boundary_test(254, 126, 0, 1, 1, 1, 16);
    do_set_velocity_on_boundary_test(254, 62, 0, 1, 1, 4, 8);
    do_set_velocity_on_boundary_test(62, 254, 0, 1, 1, 4, 8);

// --- NON-SQUARE AREA, MULTIPLE PARTITIONS, NOT FIRST LOCALITY --- //
    do_set_velocity_on_boundary_test(62, 30, 2, 2, 2, 8, 8);
    do_set_velocity_on_boundary_test(30, 62, 2, 2, 2, 8, 8);
    do_set_velocity_on_boundary_test(510, 126, 1, 1, 2, 8, 16);
    do_set_velocity_on_boundary_test(254, 510, 1, 1, 2, 8, 16);
    do_set_velocity_on_boundary_test(510, 126, 1, 2, 1, 8, 16);
    do_set_velocity_on_boundary_test(126, 510, 1, 2, 1, 8, 16);
    do_set_velocity_on_boundary_test(510, 126, 3, 2, 2, 1, 8);
    do_set_velocity_on_boundary_test(510, 254, 3, 2, 2, 1, 8);
    do_set_velocity_on_boundary_test(510, 254, 7, 4, 4, 4, 8);
    do_set_velocity_on_boundary_test(254, 510, 7, 4, 4, 4, 8);*/


    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}
