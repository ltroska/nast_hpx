#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include "grid/partition.hpp"
#include "computation/with_for_each.hpp"
#include "computation/stencils.hpp"
#include "util/cell.hpp"
#include "test_helpers.hpp"

void check_compute_rhs(scalar_grid_type const& rhs_grid, vector_grid_type const& fg_grid, index_grid_type index, computation::parameters p, RealType dt, std::string msg)
{
    using namespace computation;

    std::vector<std::vector<scalar_data> > rhs_data;

    rhs_data.resize(p.num_partitions_x);
    for (uint k = 0; k < p.num_partitions_x ; k++)
    {
        rhs_data[k].resize(p.num_partitions_y);
        for (uint l = 0; l < p.num_partitions_y ; l++)
        {
            scalar_data base = rhs_grid[l * p.num_partitions_x + k].get_data(CENTER).get();
            rhs_data[k][l] = scalar_data(base);
        }
    }

    uint global_i, global_j, size_x, size_y;

    for (uint l = 1; l < p.num_partitions_y-1; l++)
    {
        for (uint k = 1; k < p.num_partitions_x-1; k++)
        {
            size_x = rhs_data[k][l].size_x();
            size_y = rhs_data[k][l].size_y();

            for (uint j = 0; j < size_y; j++)
            {
                for (uint i = 0; i < size_x; i++)
                {
                    global_i = index[l * p.num_partitions_x + k].first + i;
                    global_j = index[l * p.num_partitions_x + k].second + j;


                    vector_data fg_center_data = fg_grid[l * p.num_partitions_x + k].get_data(CENTER).get();
                    vector_data fg_left_data = fg_grid[l * p.num_partitions_x + k-1].get_data(LEFT).get();
                    vector_data fg_bottom_data = fg_grid[(l-1) * p.num_partitions_x + k].get_data(BOTTOM).get();

                    scalar_cell const rhs_cell = rhs_data[k][l].get_cell(i, j);
                    vector_cell const fg_cell = fg_center_data.get_cell(i, j);

                    vector_cell const fg_left = get_neighbor_cell(fg_center_data, fg_left_data, fg_left_data, fg_bottom_data,
                                                                fg_bottom_data, fg_bottom_data, fg_bottom_data, fg_bottom_data, fg_bottom_data, i, j, LEFT);

                    vector_cell const fg_bottom = get_neighbor_cell(fg_center_data, fg_left_data, fg_left_data, fg_bottom_data,
                                                                fg_bottom_data, fg_bottom_data, fg_bottom_data, fg_bottom_data, fg_bottom_data, i, j, BOTTOM);


                    if (in_range(1, p.i_max, 1, p.j_max, global_i, global_j))
                    {
                        std::string ident = " partition/cell " + std::to_string(k) + " " + std::to_string(l) + " " + std::to_string(i) + " " + std::to_string(j) + " ";

                        HPX_ASSERT_MSG(rhs_cell.value == 1./dt*( (fg_cell.first - fg_left.first)/p.dx + (fg_cell.second - fg_bottom.second)/p.dy),
                                        (msg + ident + expected_string(1./dt*( (fg_cell.first - fg_left.first)/p.dx + (fg_cell.second - fg_bottom.second)/p.dy), rhs_cell.value)).c_str());
                    }
                }
            }
        }
    }
}

void do_compute_rhs_test(uint i_max, uint j_max, uint locality_id, uint localities_x, uint localities_y, uint i_res, uint j_res, RealType re = 1000, RealType dt = 1)
{

        computation::parameters params;
        params.i_max = i_max;
        params.j_max = j_max;
        params.num_cells_per_partition_x = i_res;
        params.num_cells_per_partition_y = j_res;

        params.num_partitions_x = ((i_max + 2) / localities_x) / i_res + 2;
        params.num_partitions_y = ((j_max + 2) / localities_y) / j_res + 2;
        params.re = re;
        params.dx = 0.25;
        params.dy = 0.25;
        params.alpha = 0.9;

        grid_maker maker = grid_maker(localities_x, localities_y, locality_id, params);

        index_grid_type index;
        vector_grid_type fg_grid, uv_grid;
        scalar_grid_type rhs_grid;

        maker.make_index_grid(index);
        maker.make_random_grid(uv_grid);
        maker.make_vector_grid(fg_grid, 0);
        maker.make_scalar_grid(rhs_grid, 0);


        computation::with_for_each strat(index, params);
        strat.compute_fg(fg_grid, uv_grid, dt);
        strat.compute_rhs(rhs_grid, fg_grid, dt);

        std::string msg = "\nfailed with settings " + std::to_string(i_max) + " " + std::to_string(j_max) + " " + std::to_string(locality_id) + " "
                            + std::to_string(localities_x) + " " + std::to_string(localities_y) + " " + std::to_string(i_res) + " " + std::to_string(j_res) + " ";

        check_compute_rhs(rhs_grid, fg_grid, index, params, dt, msg);
}

int hpx_main(int argc, char* argv[])
{
// --- SQUARE AREA, ONE PARTITION --- //
    do_compute_rhs_test(14, 14, 0, 1, 1, 16, 16);
    do_compute_rhs_test(62, 62, 0, 1, 1, 64, 64);

// --- SQUARE AREA, MULTIPLE PARTITIONS --- //
    do_compute_rhs_test(62, 62, 0, 1, 1, 8, 8);
    do_compute_rhs_test(64, 64, 0, 1, 1, 8, 16);
    do_compute_rhs_test(64, 64, 0, 1, 1, 32, 16);
    do_compute_rhs_test(64, 64, 0, 1, 1, 16, 32);

// --- SQUARE AREA, MULTIPLE PARTITIONS, NOT FIRST LOCALITY --- //
    do_compute_rhs_test(30, 30, 2, 2, 2, 8, 8);
    do_compute_rhs_test(128, 128, 1, 1, 2, 64, 32);
   // do_compute_rhs_test(510, 510, 1, 2, 1, 64, 32);
//    do_compute_rhs_test(1024, 1024, 3, 2, 2, 1, 8);
//    do_compute_rhs_test(1024, 1024, 7, 4, 4, 4, 8);

/*
// --- NON-SQUARE AREA, ONE PARTITION --- //
    do_compute_rhs_test(20, 40, 0, 1, 1, 1, 1);
    do_compute_rhs_test(40, 20, 0, 1, 1, 1, 1);
    do_compute_rhs_test(64, 512, 0, 1, 1, 1, 1);
    do_compute_rhs_test(512, 64, 0, 1, 1, 1, 1);

// --- NON-SQUARE AREA, MULTIPLE PARTITIONS --- //
    do_compute_rhs_test(64, 32, 0, 1, 1, 8, 8);
    do_compute_rhs_test(32, 64, 0, 1, 1, 8, 8);
    do_compute_rhs_test(128, 256, 0, 1, 1, 8, 16);
    do_compute_rhs_test(256, 128, 0, 1, 1, 8, 16);
    do_compute_rhs_test(256, 128, 0, 1, 1, 1, 16);
    do_compute_rhs_test(256, 128, 0, 1, 1, 1, 16);
    do_compute_rhs_test(256, 64, 0, 1, 1, 4, 8);
    do_compute_rhs_test(64, 256, 0, 1, 1, 4, 8);

// --- NON-SQUARE AREA, MULTIPLE PARTITIONS, NOT FIRST LOCALITY --- //
    do_compute_rhs_test(64, 32, 2, 2, 2, 8, 8);
    do_compute_rhs_test(32, 64, 2, 2, 2, 8, 8);
    do_compute_rhs_test(512, 128, 1, 1, 2, 8, 16);
    do_compute_rhs_test(256, 512, 1, 1, 2, 8, 16);
    do_compute_rhs_test(512, 128, 1, 2, 1, 8, 16);
    do_compute_rhs_test(128, 512, 1, 2, 1, 8, 16);
    do_compute_rhs_test(512, 128, 3, 2, 2, 1, 8);
    do_compute_rhs_test(512, 256, 3, 2, 2, 1, 8);
    do_compute_rhs_test(512, 256, 7, 4, 4, 4, 8);
    do_compute_rhs_test(256, 512, 7, 4, 4, 4, 8);*/


    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}

