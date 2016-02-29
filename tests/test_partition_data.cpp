#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include "grid/partition_data.hpp"
#include "util/cell.hpp"

int hpx_main(int argc, char* argv[])
{
    //test empty
    grid::partition_data<> empty_pdata;
    HPX_ASSERT_MSG(empty_pdata.size_x() == 0, "\nfailed empty data size_x");
    HPX_ASSERT_MSG(empty_pdata.size_y() == 0, "\nfailed empty data size_y");
    HPX_ASSERT_MSG(empty_pdata.size() == 0, "\nfailed empty data size");

    //test non-empty square
    grid::partition_data<> sq_pdata(7, 7, 4);
    HPX_ASSERT_MSG(sq_pdata.size_x() == 7, "\nfailed empty data size_x");
    HPX_ASSERT_MSG(sq_pdata.size_y() == 7, "\nfailed empty data size_y");
    HPX_ASSERT_MSG(sq_pdata.size() == 49, "\nfailed empty data size");
    HPX_ASSERT_MSG(sq_pdata[15] == 4, "\nfailed empty data value");

    //test non-empty non-square
    grid::partition_data<> nsq_pdata(5, 8, -1);
    HPX_ASSERT_MSG(nsq_pdata.size_x() == 5, "\nfailed empty data size_x");
    HPX_ASSERT_MSG(nsq_pdata.size_y() == 8, "\nfailed empty data size_y");
    HPX_ASSERT_MSG(nsq_pdata.size() == 40, "\nfailed empty data size");
    HPX_ASSERT_MSG(nsq_pdata[6] == -1, "\nfailed empty data value");

    //test access
    grid::partition_data<> pdata(5, 5);

    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
            pdata[j*5+i] = j*5+i;

    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
        {
            HPX_ASSERT_MSG(pdata.get_cell(i, j) == j*5+i, "\nfailed access");
            HPX_ASSERT_MSG(pdata.get_cell_ref(i, j) == j*5+i, "\nfailed access with ref");
        }

    //test directions
    grid::partition_data<> left = grid::partition_data<>(pdata, LEFT);
    HPX_ASSERT_MSG(left.size() == 5, "\nfailed LEFT size");
    HPX_ASSERT_MSG(left.size_x() == 1, "\nfailed LEFT size_x");
    HPX_ASSERT_MSG(left.size_y() == 5, "\nfailed LEFT size_y");
    HPX_ASSERT_MSG(left[1] == 9, "\nfailed LEFT value");

    grid::partition_data<> right = grid::partition_data<>(pdata, RIGHT);
    HPX_ASSERT_MSG(right.size() == 5, "\nfailed RIGHT size");
    HPX_ASSERT_MSG(right.size_x() == 1, "\nfailed RIGHT size_x");
    HPX_ASSERT_MSG(right.size_y() == 5, "\nfailed RIGHT size_y");
    HPX_ASSERT_MSG(right[2] == 10, "\nfailed RIGHT value");

    grid::partition_data<> bottom = grid::partition_data<>(pdata, BOTTOM);
    HPX_ASSERT_MSG(bottom.size() == 5, "\nfailed BOTTOM size");
    HPX_ASSERT_MSG(bottom.size_x() == 5, "\nfailed BOTTOM size_x");
    HPX_ASSERT_MSG(bottom.size_y() == 1, "\nfailed BOTTOM size_y");
    HPX_ASSERT_MSG(bottom[3] == 23, "\nfailed BOTTOM value");

    grid::partition_data<> top = grid::partition_data<>(pdata, TOP);
    HPX_ASSERT_MSG(top.size() == 5, "\nfailed TOP size");
    HPX_ASSERT_MSG(top.size_x() == 5, "\nfailed TOP size_x");
    HPX_ASSERT_MSG(top.size_y() == 1, "\nfailed TOP size_y");
    HPX_ASSERT_MSG(top[3] == 3, "\nfailed TOP value");

    grid::partition_data<> bottom_left = grid::partition_data<>(pdata, BOTTOM_LEFT);
    HPX_ASSERT_MSG(bottom_left.size() == 1, "\nfailed BOTTOM_LEFT size");
    HPX_ASSERT_MSG(bottom_left.size_x() == 1, "\nfailed BOTTOM_LEFT size_x");
    HPX_ASSERT_MSG(bottom_left.size_y() == 1, "\nfailed BOTTOM_LEFT size_y");
    HPX_ASSERT_MSG(bottom_left[0] == 24, "\nfailed BOTTOM_LEFT value");

    grid::partition_data<> bottom_right = grid::partition_data<>(pdata, BOTTOM_RIGHT);
    HPX_ASSERT_MSG(bottom_right.size() == 1, "\nfailed BOTTOM_RIGHT size");
    HPX_ASSERT_MSG(bottom_right.size_x() == 1, "\nfailed BOTTOM_RIGHT size_x");
    HPX_ASSERT_MSG(bottom_right.size_y() == 1, "\nfailed BOTTOM_RIGHT size_y");
    HPX_ASSERT_MSG(bottom_right[0] == 20, "\nfailed BOTTOM_RIGHT value");

    grid::partition_data<> top_left = grid::partition_data<>(pdata, TOP_LEFT);
    HPX_ASSERT_MSG(top_left.size() == 1, "\nfailed TOP_LEFT size");
    HPX_ASSERT_MSG(top_left.size_x() == 1, "\nfailed TOP_LEFT size_x");
    HPX_ASSERT_MSG(top_left.size_y() == 1, "\nfailed TOP_LEFT size_y");
    HPX_ASSERT_MSG(top_left[0] == 4, "\nfailed TOP_LEFT value");

    grid::partition_data<> top_right = grid::partition_data<>(pdata, TOP_RIGHT);
    HPX_ASSERT_MSG(top_right.size() == 1, "\nfailed TOP_RIGHT size");
    HPX_ASSERT_MSG(top_right.size_x() == 1, "\nfailed TOP_RIGHT size_x");
    HPX_ASSERT_MSG(top_right.size_y() == 1, "\nfailed TOP_RIGHT size_y");
    HPX_ASSERT_MSG(top_right[0] == 0, "\nfailed TOP_RIGHT value");

    grid::partition_data<> center = grid::partition_data<>(pdata, CENTER);
    HPX_ASSERT_MSG(center.size() == 25, "\nfailed CENTER size");
    HPX_ASSERT_MSG(center.size_x() == 5, "\nfailed CENTER size_x");
    HPX_ASSERT_MSG(center.size_y() == 5, "\nfailed CENTER size_y");
    HPX_ASSERT_MSG(center.get_cell(3, 3) == 18, "\nfailed CENTER value");

    //test non-empty square with custom type
    grid::partition_data<vector_cell> sq_pdata_cell(7, 7, 4);
    HPX_ASSERT_MSG(sq_pdata_cell.size_x() == 7, "\nfailed vector_cell size_x");
    HPX_ASSERT_MSG(sq_pdata_cell.size_y() == 7, "\nfailed vector_cell size_y");
    HPX_ASSERT_MSG(sq_pdata_cell.size() == 49, "\nfailed vector_cell size");
    HPX_ASSERT_MSG(sq_pdata_cell[15].first == 4 && sq_pdata_cell[15].second == 4, "\nfailed vector_cell value");

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}
