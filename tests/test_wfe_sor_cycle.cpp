#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include "grid/types.hpp"
#include "computation/with_for_each.hpp"
#include "computation/stencils.hpp"
#include "test_helpers.hpp"
#include "util/helpers.hpp"

void do_sor_cycle_test(uint i_max, uint j_max, uint locality_id, uint localities_x, uint localities_y, uint i_res, uint j_res, uint iter_max = 100, RealType re = 1000)
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
        params.omega = 1.7;

        grid_maker maker = grid_maker(localities_x, localities_y, locality_id, params);

        index_grid_type index;
        scalar_grid_type p_grid, p_grid2, rhs_grid;

        maker.make_index_grid(index);
        maker.make_random_grid(rhs_grid);
        maker.make_scalar_grid(p_grid, 0);
        maker.make_scalar_grid(p_grid2, 0);

        RealType eps = 1e-4;
        RealType new_res = 0, old_res = 1e6;

        std::string msg = "\nfailed with settings " + std::to_string(i_max) + " " + std::to_string(j_max) + " " + std::to_string(locality_id) + " "
                            + std::to_string(localities_x) + " " + std::to_string(localities_y) + " " + std::to_string(i_res) + " " + std::to_string(j_res) + " ";

        computation::with_for_each strat;

        for (uint i = 0; i < 100; i++)
        {
            for (uint k = 1; k < params.num_partitions_x -1 ; k++)
            {
                for (uint l = 1; l < params.num_partitions_y -1; l++)
                {
                    uint global_i = index[l * params.num_partitions_x + k].first;
                    uint global_j = index[l * params.num_partitions_x + k].second;

                    scalar_data base = p_grid[l * params.num_partitions_x + k].get_data(CENTER).get();

                    strat.set_pressure_on_boundary(base, global_i, global_j, params.i_max, params.j_max);

                    p_grid[l * params.num_partitions_x + k] = scalar_partition(hpx::find_here(), base);
                }
            }

            for (uint k = 1; k < params.num_partitions_x -1 ; k++)
            {
                for (uint l = 1; l < params.num_partitions_y -1; l++)
                {
                    uint global_i = index[l * params.num_partitions_x + k].first;
                    uint global_j = index[l * params.num_partitions_x + k].second;

                    scalar_data base = p_grid[l * params.num_partitions_x + k].get_data(CENTER).get();
                    scalar_data rhs_center = rhs_grid[l * params.num_partitions_x + k].get_data(CENTER).get();
                    scalar_data p_left = p_grid[l * params.num_partitions_x + k -1].get_data(LEFT).get();
                    scalar_data p_right = p_grid[l * params.num_partitions_x + k +1].get_data(RIGHT).get();
                    scalar_data p_bottom = p_grid[(l-1) * params.num_partitions_x + k].get_data(BOTTOM).get();
                    scalar_data p_top = p_grid[(l+1) * params.num_partitions_x + k].get_data(TOP).get();

                    strat.sor_cycle(base, p_left, p_right, p_bottom, p_top, rhs_center,
                                    index[l * params.num_partitions_x + k].first,
                                    index[l * params.num_partitions_x + k].second, params.i_max, params.j_max,
                                    params.omega, params.dx, params.dy);

                    p_grid2[l * params.num_partitions_x + k] = scalar_partition(hpx::find_here(), base);
                }
            }

            for (uint k = 1; k < params.num_partitions_x -1 ; k++)
            {
                for (uint l = 1; l < params.num_partitions_y -1; l++)
                {
                    uint global_i = index[l * params.num_partitions_x + k].first;
                    uint global_j = index[l * params.num_partitions_x + k].second;

                    scalar_data base = p_grid2[l * params.num_partitions_x + k].get_data(CENTER).get();
                    scalar_data rhs_center = rhs_grid[l * params.num_partitions_x + k].get_data(CENTER).get();
                    scalar_data p_left = p_grid2[l * params.num_partitions_x + k -1].get_data(LEFT).get();
                    scalar_data p_right = p_grid2[l * params.num_partitions_x + k +1].get_data(RIGHT).get();
                    scalar_data p_bottom = p_grid2[(l-1) * params.num_partitions_x + k].get_data(BOTTOM).get();
                    scalar_data p_top = p_grid2[(l+1) * params.num_partitions_x + k].get_data(TOP).get();

                    new_res += strat.compute_residual(base, p_left, p_right, p_bottom, p_top, rhs_center,
                                    index[l * params.num_partitions_x + k].first,
                                    index[l * params.num_partitions_x + k].second, params.i_max, params.j_max,
                                    params.dx, params.dy);
                }
            }

            HPX_ASSERT_MSG(new_res - old_res <= eps, (msg + " old res: " + std::to_string(old_res) + " new res: " + std::to_string(new_res)).c_str());
            old_res = new_res;
            new_res = 0;
        }
}

int hpx_main(int argc, char* argv[])
{
    for (uint i = 0; i < 10; i++)
    do_sor_cycle_test(126, 126, 0, 1, 1, 128, 128);

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}


