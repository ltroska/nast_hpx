#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include "grid/partition.hpp"
#include "computation/custom_grain_size.hpp"
#include "util/cell.hpp"
#include "test_helpers.hpp"

void check_update_velocities(std::pair<RealType, RealType> max_uv, vector_grid_type const& uv_grid, vector_grid_type const& fg_grid, scalar_grid_type const& p_grid,
                                index_grid_type index, computation::parameters p, RealType dt, std::string msg)
{
    using namespace computation;

    std::vector<std::vector<vector_data> > uv_data;

    uv_data.resize(p.num_partitions_x);
    for (uint k = 0; k < p.num_partitions_x ; k++)
    {
        uv_data[k].resize(p.num_partitions_y);
        for (uint l = 0; l < p.num_partitions_y ; l++)
        {
            vector_data base = uv_grid[l * p.num_partitions_x + k].get_data(CENTER).get();
            uv_data[k][l] = vector_data(base);
        }
    }

    std::vector<std::vector<vector_data> > fg_data;

    fg_data.resize(p.num_partitions_x);
    for (uint k = 0; k < p.num_partitions_x ; k++)
    {
        fg_data[k].resize(p.num_partitions_y);
        for (uint l = 0; l < p.num_partitions_y ; l++)
        {
            vector_data base = fg_grid[l * p.num_partitions_x + k].get_data(CENTER).get();
            fg_data[k][l] = vector_data(base);
        }
    }

    uint global_i, global_j, size_x, size_y;

    for (uint l = 1; l < p.num_partitions_y-1; l++)
    {
        for (uint k = 1; k < p.num_partitions_x-1; k++)
        {
            size_x = fg_data[k][l].size_x();
            size_y = fg_data[k][l].size_y();

            for (uint j = 0; j < size_y; j++)
            {
                for (uint i = 0; i < size_x; i++)
                {
                    scalar_data const p_center = p_grid[l*p.num_partitions_x + k].get_data(CENTER).get();
                    scalar_data const p_right_data = p_grid[l*p.num_partitions_x + k+1].get_data(RIGHT).get();
                    scalar_data const p_top_data = p_grid[(l+1)*p.num_partitions_x + k].get_data(TOP).get();

                    global_i = index[l * p.num_partitions_x + k].first + i;
                    global_j = index[l * p.num_partitions_x + k].second + j;

                    vector_cell const uv_cell = uv_data[k][l].get_cell(i, j);
                    vector_cell const fg_cell = fg_data[k][l].get_cell(i, j);
                    scalar_cell const p_cell = p_center.get_cell(i, j);
                    scalar_cell const p_right = get_neighbor_cell(p_center, p_right_data, p_right_data, p_right_data,
                                                                    p_top_data, p_top_data, p_top_data, p_top_data, p_top_data, i, j, RIGHT);
                    scalar_cell const p_top = get_neighbor_cell(p_center, p_right_data, p_right_data, p_right_data,
                                                                    p_top_data, p_top_data, p_top_data, p_top_data, p_top_data, i, j, TOP);

                    std::string ident = " partition/cell " + std::to_string(k) + " " + std::to_string(l) + " " + std::to_string(i) + " " + std::to_string(j) + " ";

                    if (in_range(1, p.i_max - 1, 1, p.j_max, global_i, global_j))
                    {
                        RealType exp = fg_cell.first - dt/p.dx*(p_right.value - p_cell.value);
                        HPX_ASSERT_MSG(exp == uv_cell.first, (msg + ident + expected_string(exp, uv_cell.first)).c_str());
                    }

                    if (in_range(1, p.i_max, 1, p.j_max - 1, global_i, global_j))
                    {
                        RealType exp = fg_cell.second - dt/p.dy*(p_top.value - p_cell.value);
                        HPX_ASSERT_MSG(exp == uv_cell.second, (msg + ident + expected_string(exp, uv_cell.second)).c_str());
                    }

                    HPX_ASSERT_MSG(uv_cell.first <= max_uv.first, (msg + ident + " not max u " + expected_string(uv_cell.first, max_uv.first)).c_str());
                    HPX_ASSERT_MSG(uv_cell.second <= max_uv.second, (msg + ident + " not max v " + expected_string(uv_cell.second, max_uv.second)).c_str());
                }
            }
        }
    }
}


void do_update_velocities_test(uint i_max, uint j_max, uint locality_id, uint localities_x, uint localities_y, uint i_res,
                                    uint j_res, RealType re = 1000, RealType dt = 1)
{
        uint num_partitions_x = i_res + 2;
        uint num_partitions_y = j_res + 2;

        computation::parameters params;
        params.num_partitions_x = num_partitions_x;
        params.num_partitions_y = num_partitions_y;
        params.i_max = i_max;
        params.j_max = j_max;
        params.num_cells_per_partition_x = ((i_max + 2) / localities_x) / i_res;
        params.num_cells_per_partition_y = ((j_max + 2) / localities_y) / j_res;
        params.re = re;
        params.dx = 0.25;
        params.dy = 0.25;
        params.alpha = 0.9;

        grid_maker maker = grid_maker(localities_x, localities_y, locality_id, params);

        index_grid_type index;
        vector_grid_type fg_grid, uv_grid;
        scalar_grid_type p_grid;

        maker.make_index_grid(index);
        maker.make_random_grid(p_grid);
        maker.make_random_grid(fg_grid);
        maker.make_vector_grid(uv_grid, 0);

        computation::custom_grain_size strat(index, params);
        std::pair<RealType, RealType> max_uv = strat.update_velocities(uv_grid, fg_grid, p_grid, dt).get();

        std::string msg = "\nfailed with settings " + std::to_string(i_max) + " " + std::to_string(j_max) + " " + std::to_string(locality_id) + " "
                            + std::to_string(localities_x) + " " + std::to_string(localities_y) + " " + std::to_string(i_res) + " " + std::to_string(j_res) + " ";

        check_update_velocities(max_uv, uv_grid, fg_grid, p_grid, index, params, dt, msg);
}

int hpx_main(int argc, char* argv[])
{
// --- SQUARE AREA, ONE PARTITION --- //
    do_update_velocities_test(6, 6, 0, 1, 1, 1, 1);
    do_update_velocities_test(62, 62, 0, 1, 1, 1, 1);

// --- SQUARE AREA, MULTIPLE PARTITIONS --- //
    do_update_velocities_test(62, 62, 0, 1, 1, 8, 8);
    do_update_velocities_test(254, 254, 0, 1, 1, 8, 16);
    do_update_velocities_test(254, 254, 0, 1, 1, 1, 8);
    do_update_velocities_test(254, 254, 0, 1, 1, 4, 8);

// --- SQUARE AREA, MULTIPLE PARTITIONS, NOT FIRST LOCALITY --- //
    do_update_velocities_test(30, 30, 2, 2, 2, 8, 8);
    do_update_velocities_test(510, 510, 1, 1, 2, 8, 16);
    do_update_velocities_test(510, 510, 1, 2, 1, 8, 16);
//    do_update_velocities_test(1022, 1022, 3, 2, 2, 1, 8);
//    do_update_velocities_test(1022, 1022, 7, 4, 4, 4, 8);

/*
// --- NON-SQUARE AREA, ONE PARTITION --- //
    do_update_velocities_test(14, 30, 0, 1, 1, 1, 1);
    do_update_velocities_test(30, 14, 0, 1, 1, 1, 1);
    do_update_velocities_test(62, 510, 0, 1, 1, 1, 1);
    do_update_velocities_test(510, 62, 0, 1, 1, 1, 1);

// --- NON-SQUARE AREA, MULTIPLE PARTITIONS --- //
    do_update_velocities_test(62, 30, 0, 1, 1, 8, 8);
    do_update_velocities_test(30, 62, 0, 1, 1, 8, 8);
    do_update_velocities_test(126, 254, 0, 1, 1, 8, 16);
    do_update_velocities_test(254, 126, 0, 1, 1, 8, 16);
    do_update_velocities_test(254, 126, 0, 1, 1, 1, 16);
    do_update_velocities_test(254, 126, 0, 1, 1, 1, 16);
    do_update_velocities_test(254, 62, 0, 1, 1, 4, 8);
    do_update_velocities_test(62, 254, 0, 1, 1, 4, 8);

// --- NON-SQUARE AREA, MULTIPLE PARTITIONS, NOT FIRST LOCALITY --- //
    do_update_velocities_test(62, 30, 2, 2, 2, 8, 8);
    do_update_velocities_test(30, 62, 2, 2, 2, 8, 8);
    do_update_velocities_test(510, 126, 1, 1, 2, 8, 16);
    do_update_velocities_test(254, 510, 1, 1, 2, 8, 16);
    do_update_velocities_test(510, 126, 1, 2, 1, 8, 16);
    do_update_velocities_test(126, 510, 1, 2, 1, 8, 16);
    do_update_velocities_test(510, 126, 3, 2, 2, 1, 8);
    do_update_velocities_test(510, 254, 3, 2, 2, 1, 8);
    do_update_velocities_test(510, 254, 7, 4, 4, 4, 8);
    do_update_velocities_test(254, 510, 7, 4, 4, 4, 8);*/

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}
