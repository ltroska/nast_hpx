#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include "grid/partition_data.hpp"
#include "util/cell.hpp"

int hpx_main(int argc, char* argv[])
{
    //test empty
    grid::partition_data<> empty_pdata;
    HPX_ASSERT(empty_pdata.size_x() == 0);
    HPX_ASSERT(empty_pdata.size_y() == 0);
    HPX_ASSERT(empty_pdata.size() == 0);

    //test non-empty square
    grid::partition_data<> sq_pdata(7, 7, 4);
    HPX_ASSERT(sq_pdata.size_x() == 7);
    HPX_ASSERT(sq_pdata.size_y() == 7);
    HPX_ASSERT(sq_pdata.size() == 49);
    HPX_ASSERT(sq_pdata[15] == 4);

    //test non-empty non-square
    grid::partition_data<> nsq_pdata(5, 8, -1);
    HPX_ASSERT(nsq_pdata.size_x() == 5);
    HPX_ASSERT(nsq_pdata.size_y() == 8);
    HPX_ASSERT(nsq_pdata.size() == 40);
    HPX_ASSERT(nsq_pdata[6] == -1);

    //test access
    grid::partition_data<> pdata(5, 5);

    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
            pdata[j*5+i] = j*5+i;

    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
        {
            HPX_ASSERT(pdata.get_cell(i, j) == j*5+i);
            HPX_ASSERT(pdata.get_cell_ref(i, j) == j*5+i);
        }

    //test directions
    grid::partition_data<> left = grid::partition_data<>(pdata, LEFT);
    HPX_ASSERT(left.size() == 5);
    HPX_ASSERT(left.size_x() == 1);
    HPX_ASSERT(left.size_y() == 5);
    HPX_ASSERT(left[1] == 9);

    grid::partition_data<> right = grid::partition_data<>(pdata, RIGHT);
    HPX_ASSERT(right.size() == 5);
    HPX_ASSERT(right.size_x() == 1);
    HPX_ASSERT(right.size_y() == 5);
    HPX_ASSERT(right[2] == 10);

    grid::partition_data<> bottom = grid::partition_data<>(pdata, BOTTOM);
    HPX_ASSERT(bottom.size() == 5);
    HPX_ASSERT(bottom.size_x() == 5);
    HPX_ASSERT(bottom.size_y() == 1);
    HPX_ASSERT(bottom[3] == 23);

    grid::partition_data<> top = grid::partition_data<>(pdata, TOP);
    HPX_ASSERT(top.size() == 5);
    HPX_ASSERT(top.size_x() == 5);
    HPX_ASSERT(top.size_y() == 1);
    HPX_ASSERT(top[3] == 3);

    grid::partition_data<> bottom_left = grid::partition_data<>(pdata, BOTTOM_LEFT);
    HPX_ASSERT(bottom_left.size() == 1);
    HPX_ASSERT(bottom_left.size_x() == 1);
    HPX_ASSERT(bottom_left.size_y() == 1);
    HPX_ASSERT(bottom_left[0] == 24);

    grid::partition_data<> bottom_right = grid::partition_data<>(pdata, BOTTOM_RIGHT);
    HPX_ASSERT(bottom_right.size() == 1);
    HPX_ASSERT(bottom_right.size_x() == 1);
    HPX_ASSERT(bottom_right.size_y() == 1);
    HPX_ASSERT(bottom_right[0] == 20);

    grid::partition_data<> top_left = grid::partition_data<>(pdata, TOP_LEFT);
    HPX_ASSERT(top_left.size() == 1);
    HPX_ASSERT(top_left.size_x() == 1);
    HPX_ASSERT(top_left.size_y() == 1);
    HPX_ASSERT(top_left[0] == 4);

    grid::partition_data<> top_right = grid::partition_data<>(pdata, TOP_RIGHT);
    HPX_ASSERT(top_right.size() == 1);
    HPX_ASSERT(top_right.size_x() == 1);
    HPX_ASSERT(top_right.size_y() == 1);
    HPX_ASSERT(top_right[0] == 0);

    grid::partition_data<> center = grid::partition_data<>(pdata, CENTER);
    HPX_ASSERT(center.size() == 25);
    HPX_ASSERT(center.size_x() == 5);
    HPX_ASSERT(center.size_y() == 5);
    HPX_ASSERT(center.get_cell(3, 3) == 18);

    //test non-empty square with custom type
    grid::partition_data<vector_cell> sq_pdata_cell(7, 7, 4);
    HPX_ASSERT(sq_pdata_cell.size_x() == 7);
    HPX_ASSERT(sq_pdata_cell.size_y() == 7);
    HPX_ASSERT(sq_pdata_cell.size() == 49);
    HPX_ASSERT(sq_pdata_cell[15].first == 4 && sq_pdata_cell[15].second == 4);

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}
