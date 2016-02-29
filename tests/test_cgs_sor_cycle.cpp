#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include "grid/partition.hpp"
#include "computation/custom_grain_size.hpp"
#include "computation/stencils.hpp"
#include "util/cell.hpp"
#include "test_helpers.hpp"

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
        scalar_grid_type p_grid, rhs_grid;

        maker.make_index_grid(index);
        maker.make_random_grid(rhs_grid);
        maker.make_scalar_grid(p_grid, 0);

        RealType eps = 1e-4;
        RealType new_res, old_res = 1e6;

        std::string msg = "\nfailed with settings " + std::to_string(i_max) + " " + std::to_string(j_max) + " " + std::to_string(locality_id) + " "
                            + std::to_string(localities_x) + " " + std::to_string(localities_y) + " " + std::to_string(i_res) + " " + std::to_string(j_res) + " ";

        computation::custom_grain_size strat(index, params);
        for (uint i = 0; i < 100; i++)
        {
            strat.set_pressure_on_boundary(p_grid);
            strat.sor_cycle(p_grid, rhs_grid);
            new_res = strat.compute_residual(p_grid, rhs_grid).get();
            HPX_ASSERT_MSG(new_res - old_res <= eps, (msg + " old res: " + std::to_string(old_res) + " new res: " + std::to_string(new_res)).c_str());
            old_res = new_res;
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


