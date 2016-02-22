#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include "grid/partition.hpp"
#include "computation/custom_grain_size.hpp"
#include "computation/stencils.hpp"
#include "util/cell.hpp"
#include "test_helpers.hpp"

void check_compute_fg(vector_grid_type const& fg_grid, vector_grid_type const& uv_grid, index_grid_type index, computation::parameters p, RealType dt, std::string msg)
{
    using namespace computation;

    std::vector<std::vector<vector_data > > fg_data;

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

    std::vector<std::vector<vector_data > > uv_data;

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
                    global_i = index[l * p.num_partitions_x + k].first + i;
                    global_j = index[l * p.num_partitions_x + k].second + j;

                    vector_cell fg_cell = fg_data[k][l].get_cell(i, j);
                    vector_cell uv_cell = uv_data[k][l].get_cell(i, j);

                    vector_cell uv_left = get_neighbor_cell(uv_data[k][l], uv_data[k-1][l], uv_data[k+1][l], uv_data[k][l-1], uv_data[k][l+1],
                                                uv_data[k-1][l-1], uv_data[k+1][l-1], uv_data[k-1][l+1], uv_data[k+1][l+1], i, j, LEFT);
                    vector_cell uv_right = get_neighbor_cell(uv_data[k][l], uv_data[k-1][l], uv_data[k+1][l], uv_data[k][l-1], uv_data[k][l+1],
                                                uv_data[k-1][l-1], uv_data[k+1][l-1], uv_data[k-1][l+1], uv_data[k+1][l+1], i, j, RIGHT);
                    vector_cell uv_bottom = get_neighbor_cell(uv_data[k][l], uv_data[k-1][l], uv_data[k+1][l], uv_data[k][l-1], uv_data[k][l+1],
                                                uv_data[k-1][l-1], uv_data[k+1][l-1], uv_data[k-1][l+1], uv_data[k+1][l+1], i, j, BOTTOM);
                    vector_cell uv_top = get_neighbor_cell(uv_data[k][l], uv_data[k-1][l], uv_data[k+1][l], uv_data[k][l-1], uv_data[k][l+1],
                                                uv_data[k-1][l-1], uv_data[k+1][l-1], uv_data[k-1][l+1], uv_data[k+1][l+1], i, j, TOP);
                    vector_cell uv_bottomleft = get_neighbor_cell(uv_data[k][l], uv_data[k-1][l], uv_data[k+1][l], uv_data[k][l-1], uv_data[k][l+1],
                                                uv_data[k-1][l-1], uv_data[k+1][l-1], uv_data[k-1][l+1], uv_data[k+1][l+1], i, j, BOTTOM_LEFT);
                    vector_cell uv_bottomright = get_neighbor_cell(uv_data[k][l], uv_data[k-1][l], uv_data[k+1][l], uv_data[k][l-1], uv_data[k][l+1],
                                                uv_data[k-1][l-1], uv_data[k+1][l-1], uv_data[k-1][l+1], uv_data[k+1][l+1], i, j, BOTTOM_RIGHT);
                    vector_cell uv_topleft = get_neighbor_cell(uv_data[k][l], uv_data[k-1][l], uv_data[k+1][l], uv_data[k][l-1], uv_data[k][l+1],
                                                uv_data[k-1][l-1], uv_data[k+1][l-1], uv_data[k-1][l+1], uv_data[k+1][l+1], i, j, TOP_LEFT);
                    vector_cell uv_topright = get_neighbor_cell(uv_data[k][l], uv_data[k-1][l], uv_data[k+1][l], uv_data[k][l-1], uv_data[k][l+1],
                                                uv_data[k-1][l-1], uv_data[k+1][l-1], uv_data[k-1][l+1], uv_data[k+1][l+1], i, j, TOP_RIGHT);

                    std::string ident = std::to_string(k) + " " + std::to_string(l) + " " + std::to_string(i) + " " + std::to_string(j);
                    std::cout << fg_cell ;

                    if (in_range(1, p.i_max-1, 1, p.j_max, global_i, global_j))
                    {
                        RealType res = uv_cell.first
                                        + dt* (1./p.re * (second_derivative_fwd_bkwd_x(uv_right.first, uv_cell.first, uv_left.first, p.dx)
                                                            + second_derivative_fwd_bkwd_y(uv_top.first, uv_cell.first, uv_bottom.first, p.dy))

                                                - first_derivative_of_square_x(uv_right.first, uv_cell.first, uv_left.first, p.dx, p.alpha)
                                                - first_derivative_of_product_y(uv_right.second, uv_cell.second, uv_bottom.second, uv_bottomright.second,
                                                                                uv_bottom.first, uv_cell.first, uv_top.first, p.dy, p.alpha)
                                                    );

                        std::cout  << " " <<res ;


                        HPX_ASSERT_MSG(res == fg_cell.first,
                                        (msg + ident + expected_string(res, fg_cell.first)).c_str());
                    }

                    if (in_range(1, p.i_max, 1, p.j_max-1, global_i, global_j))
                    {
                        RealType res = uv_cell.second
                                        + dt*(1./p.re * (second_derivative_fwd_bkwd_x(uv_right.second,  uv_cell.second, uv_left.second, p.dx)
                                                            + second_derivative_fwd_bkwd_y(uv_top.second, uv_cell.second, uv_bottom.second, p.dy))
                                                    - first_derivative_of_product_x(uv_left.first, uv_cell.first, uv_top.first, uv_topleft.first,
                                                        uv_left.second, uv_cell.second, uv_right.second, p.dx, p.alpha)
                                                    );
                        std::cout  << " " <<res ;

                        HPX_ASSERT_MSG(res == fg_cell.second, (msg + ident + expected_string(res, fg_cell.second)).c_str());
                    }
                    std::cout << std::endl;
                }
            }
        }
    }
}

void do_compute_fg_test(uint i_max, uint j_max, uint locality_id, uint localities_x, uint localities_y, uint i_res, uint j_res, RealType re = 1000, RealType dt = 1)
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
        params.re = re;
        params.dx = 0.25;
        params.dy = 0.25;
        params.alpha = 0.9;

        grid_maker maker = grid_maker(localities_x, localities_y, locality_id, params);

        index_grid_type index;
        vector_grid_type fg_grid, uv_grid;

        maker.make_index_grid(index);
        maker.make_vector_grid(fg_grid, 0);
        maker.make_vector_grid(uv_grid, 2);

       print_grid(fg_grid, params);
        print_grid(uv_grid, params);

        computation::custom_grain_size strat(index, params);
        strat.compute_fg(fg_grid, uv_grid, dt);
        print_grid(fg_grid, params);

        std::string msg = "\nfailed with settings " + std::to_string(i_max) + " " + std::to_string(j_max) + " " + std::to_string(locality_id) + " "
                            + std::to_string(localities_x) + " " + std::to_string(localities_y) + " " + std::to_string(i_res) + " " + std::to_string(j_res) + " ";

        check_compute_fg(fg_grid, uv_grid, index, params, dt, msg);
}

int hpx_main(int argc, char* argv[])
{
    do_compute_fg_test(4, 4, 0, 1, 1, 1, 1);

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}
