
#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include "grid/types.hpp"
#include "util/helpers.hpp"
#include "test_helpers.hpp"

void do_get_neighbor_cell_test(uint k, uint l, uint i, uint j, uint i_max, uint j_max, uint i_res, uint j_res,
                                    uint locality_id = 0, uint localities_x = 1, uint localities_y = 1)
{
    uint num_partitions_x = i_res + 2;
    uint num_partitions_y = j_res + 2;

    k++;
    l++;

    uint num_cells_per_partition_x = ((i_max + 2) / localities_x) / i_res;
    uint num_cells_per_partition_y = ((j_max + 2) / localities_y) / j_res;

    vector_grid_type grid;
    grid_maker maker = grid_maker(localities_x, localities_y, locality_id,
                                    num_cells_per_partition_x,
                                    num_cells_per_partition_y,
                                    num_partitions_x, num_partitions_y);
    maker.make_neighbor_test_grid(grid);
    
    vector_data center_data = grid[(l) * num_partitions_y + k].get_data(CENTER).get();
    vector_data left_data = grid[(l) * num_partitions_y + k - 1].get_data(LEFT).get();
    vector_data right_data = grid[(l) * num_partitions_y + k + 1].get_data(RIGHT).get();
    vector_data bottom_data = grid[(l-1) * num_partitions_y + k].get_data(BOTTOM).get();
    vector_data top_data = grid[(l+1) * num_partitions_y + k].get_data(TOP).get();
    vector_data bottomleft_data = grid[(l-1) * num_partitions_y + k-1].get_data(BOTTOM_LEFT).get();
    vector_data bottomright_data = grid[(l-1) * num_partitions_y + k+1].get_data(BOTTOM_RIGHT).get();
    vector_data topleft_data = grid[(l+1) * num_partitions_y + k-1].get_data(TOP_LEFT).get();
    vector_data topright_data = grid[(l+1) * num_partitions_y + k+1].get_data(TOP_RIGHT).get();

    vector_cell& cell = center_data.get_cell_ref(i, j);
    vector_cell left_neighbor = get_left_neighbor(center_data, left_data, i, j);
    vector_cell right_neighbor =
        get_right_neighbor(center_data, right_data, i, j);
    vector_cell bottom_neighbor =
        get_bottom_neighbor(center_data, bottom_data, i, j);
    vector_cell top_neighbor = get_top_neighbor(center_data, top_data, i, j);
    vector_cell bottom_right_neighbor =
        get_bottomright_neighbor(center_data, bottom_data, right_data,
            bottomright_data, i, j);
    vector_cell top_left_neighbor =
        get_left_neighbor(center_data, left_data, i, j);

    std::string msg = std::to_string(k) + " " + std::to_string(l) + " " + std::to_string(i) + " " + std::to_string(j);
    HPX_ASSERT_MSG(left_neighbor.first == cell.first - 1 && left_neighbor.second == cell.second,
                    ("left_neighbor failed for " + msg + expected_string(cell, left_neighbor)).c_str() );

    HPX_ASSERT_MSG(right_neighbor.first == cell.first + 1 && right_neighbor.second == cell.second,
                    ("right_neighbor failed for " + msg + expected_string(cell, right_neighbor)).c_str() );

    HPX_ASSERT_MSG(bottom_neighbor.first == cell.first && bottom_neighbor.second == cell.second - 1,
                    ("bottom_neighbor failed for " + msg + expected_string(cell, bottom_neighbor)).c_str() );

    HPX_ASSERT_MSG(top_neighbor.first == cell.first && top_neighbor.second == cell.second + 1,
                    ("top_neighbor failed for " + msg + expected_string(cell, top_neighbor)).c_str() );

    HPX_ASSERT_MSG(bottom_right_neighbor.first == cell.first + 1 && bottom_right_neighbor.second == cell.second - 1,
                    ("bottom_right_neighbor failed for " + msg + expected_string(cell, bottom_right_neighbor)).c_str() );

    HPX_ASSERT_MSG(top_left_neighbor.first == cell.first - 1 && top_left_neighbor.second == cell.second + 1,
                    ("top_left_neighbor failed for " + msg + expected_string(cell, top_left_neighbor)).c_str() );

}


int hpx_main(int argc, char* argv[])
{
    uint i_max = 6;
    uint j_max = 6;
    uint i_res = 1;
    uint j_res = 1;

    for (uint i = 0; i < i_max; i++)
        for (uint j = 0; j < j_max; j++)
            do_get_neighbor_cell_test(0, 0, i, j, i_max, j_max, i_res, j_res);



    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}
