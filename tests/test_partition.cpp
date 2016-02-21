#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include "grid/partition.hpp"
#include "util/cell.hpp"

HPX_REGISTER_PARTITION_SERVER(RealType);

int hpx_main(int argc, char* argv[])
{
    std::size_t x = 5;
    std::size_t y = 5;
    std::size_t size = x*y;

    // ------------------------------- RealType ------------------------------- //
    {
        typedef RealType type;

        //test empty
        grid::partition<type> empty_part(hpx::find_here(), x, y);
        grid::partition_data<type> empty_pdata = empty_part.get_data(CENTER).get();

        for (int i = 0; i < size; i++)
            HPX_ASSERT_MSG(empty_pdata[i] == 0, "\nRealType, failed empty");

        //test initial value
        grid::partition<type> initial_part(hpx::find_here(), x, y, 4);
        grid::partition_data<type> initial_pdata = initial_part.get_data(CENTER).get();

        for (int i = 0; i < size; i++)
            HPX_ASSERT_MSG(initial_pdata[i] == 4, "\nRealType, failed setting initial data");

        //test base data
        grid::partition_data<type> base_pdata(5, 5);
        for (int i = 0; i < size; i++)
                base_pdata[i] = i;

        grid::partition<type> base_part(hpx::find_here(), base_pdata);
        grid::partition_data<type> base_recv_pdata = base_part.get_data(CENTER).get();
        for (int i = 0; i < size; i++)
                HPX_ASSERT_MSG(base_recv_pdata[i] == base_pdata[i], "\nRealType, failed receiving center data");

        //test directional
        grid::partition_data<type> left = base_part.get_data(LEFT).get();
        HPX_ASSERT_MSG(left.size() == y, "\nRealType, failed LEFT size");
        HPX_ASSERT_MSG(left.size_x() == 1, "\nRealType, failed LEFT size_x");
        HPX_ASSERT_MSG(left.size_y() == y, "\nRealType, failed LEFT size_y");
        HPX_ASSERT_MSG(left[1] == 9, "\nRealType, failed LEFT value");

        grid::partition_data<type> right = base_part.get_data(RIGHT).get();
        HPX_ASSERT_MSG(right.size() == y, "\nRealType, failed RIGHT size");
        HPX_ASSERT_MSG(right.size_x() == 1, "\nRealType, failed RIGHT size_x");
        HPX_ASSERT_MSG(right.size_y() == y, "\nRealType, failed RIGHT size_y");
        HPX_ASSERT_MSG(right[2] == 10, "\nRealType, failed RIGHT value");

        grid::partition_data<type> bottom = base_part.get_data(BOTTOM).get();
        HPX_ASSERT_MSG(bottom.size() == x, "\nRealType, failed BOTTOM size");
        HPX_ASSERT_MSG(bottom.size_x() == x, "\nRealType, failed BOTTOM size_x");
        HPX_ASSERT_MSG(bottom.size_y() == 1, "\nRealType, failed BOTTOM size_y");
        HPX_ASSERT_MSG(bottom[3] == 23, "\nRealType, failed BOTTOM value");

        grid::partition_data<type> top = base_part.get_data(TOP).get();
        HPX_ASSERT_MSG(top.size() == x, "\nRealType, failed TOP size");
        HPX_ASSERT_MSG(top.size_x() == x, "\nRealType, failed TOP size_x");
        HPX_ASSERT_MSG(top.size_y() == 1, "\nRealType, failed TOP size_y");
        HPX_ASSERT_MSG(top[3] == 3, "\nRealType, failed TOP value");

        grid::partition_data<type> bottom_left = base_part.get_data(BOTTOM_LEFT).get();
        HPX_ASSERT_MSG(bottom_left.size() == 1, "\nRealType, failed BOTTOM_LEFT size");
        HPX_ASSERT_MSG(bottom_left.size_x() == 1, "\nRealType, failed BOTTOM_LEFT size_x");
        HPX_ASSERT_MSG(bottom_left.size_y() == 1, "\nRealType, failed BOTTOM_LEFT size_y");
        HPX_ASSERT_MSG(bottom_left[0] == 24, "\nRealType, failed BOTTOM_LEFT value");

        grid::partition_data<type> bottom_right = base_part.get_data(BOTTOM_RIGHT).get();
        HPX_ASSERT_MSG(bottom_right.size() == 1, "\nRealType, failed BOTTOM_RIGHT size");
        HPX_ASSERT_MSG(bottom_right.size_x() == 1, "\nRealType, failed BOTTOM_RIGHT size_x");
        HPX_ASSERT_MSG(bottom_right.size_y() == 1, "\nRealType, failed BOTTOM_RIGHT size_y");
        HPX_ASSERT_MSG(bottom_right[0] == 20, "\nRealType, failed BOTTOM_RIGHT value");

        grid::partition_data<type> top_left = base_part.get_data(TOP_LEFT).get();
        HPX_ASSERT_MSG(top_left.size() == 1, "\nRealType, failed TOP_LEFT size");
        HPX_ASSERT_MSG(top_left.size_x() == 1, "\nRealType, failed TOP_LEFT size_x");
        HPX_ASSERT_MSG(top_left.size_y() == 1, "\nRealType, failed TOP_LEFT size_y");
        HPX_ASSERT_MSG(top_left[0] == 4, "\nRealType, failed TOP_LEFT value");

        grid::partition_data<type> top_right = base_part.get_data(TOP_RIGHT).get();
        HPX_ASSERT_MSG(top_right.size() == 1, "\nRealType, failed TOP_RIGHT size");
        HPX_ASSERT_MSG(top_right.size_x() == 1, "\nRealType, failed TOP_RIGHT size");
        HPX_ASSERT_MSG(top_right.size_y() == 1, "\nRealType, failed TOP_RIGHT size_y");
        HPX_ASSERT_MSG(top_right[0] == 0, "\nRealType, failed TOP_RIGHT value");
    }

    // ------------------------------- CUSTOM TYPE ------------------------------- //
    x = 10;
    y = 10;
    size = x * y;

    {
        typedef scalar_cell type;

        //test empty
        grid::partition<type> empty_part(hpx::find_here(), x, y);
        grid::partition_data<type> empty_pdata = empty_part.get_data(CENTER).get();

        for (int i = 0; i < size; i++)
            HPX_ASSERT_MSG(empty_pdata[i].value == 0, "\nscalar_cell, failed empty");

        //test initial value
        grid::partition<type> initial_part(hpx::find_here(), x, y, 4);
        grid::partition_data<type> initial_pdata = initial_part.get_data(CENTER).get();

        for (int i = 0; i < size; i++)
            HPX_ASSERT_MSG(initial_pdata[i].value == 4, "\nscalar_cell, failed setting initial data");

        //test base data
        grid::partition_data<type> base_pdata(x, y);
        for (int i = 0; i < size; i++)
                base_pdata[i].value = i;

        grid::partition<type> base_part(hpx::find_here(), base_pdata);
        grid::partition_data<type> base_recv_pdata = base_part.get_data(CENTER).get();
        for (int i = 0; i < size; i++)
                HPX_ASSERT_MSG(base_recv_pdata[i].value == base_pdata[i].value, "\nscalar_cell, failed receiving center data");

        //test directional
        grid::partition_data<type> left = base_part.get_data(LEFT).get();
        HPX_ASSERT_MSG(left.size() == y, "\nscalar_cell, failed LEFT size");
        HPX_ASSERT_MSG(left.size_x() == 1, "\nscalar_cell, failed LEFT size_x");
        HPX_ASSERT_MSG(left.size_y() == y, "\nscalar_cell, failed LEFT size_y");
        HPX_ASSERT_MSG(left[3].value == 39, "\nscalar_cell, failed LEFT value");

        grid::partition_data<type> right = base_part.get_data(RIGHT).get();
        HPX_ASSERT_MSG(right.size() == y, "\nscalar_cell, failed RIGHT size");
        HPX_ASSERT_MSG(right.size_x() == 1, "\nscalar_cell, failed RIGHT size_x");
        HPX_ASSERT_MSG(right.size_y() == y, "\nscalar_cell, failed RIGHT size_y");
        HPX_ASSERT_MSG(right[4].value == 40, "\nscalar_cell, failed RIGHT value");

        grid::partition_data<type> bottom = base_part.get_data(BOTTOM).get();
        HPX_ASSERT_MSG(bottom.size() == x, "\nscalar_cell, failed BOTTOM size");
        HPX_ASSERT_MSG(bottom.size_x() == x, "\nscalar_cell, failed BOTTOM size_x");
        HPX_ASSERT_MSG(bottom.size_y() == 1, "\nscalar_cell, failed BOTTOM size_y");
        HPX_ASSERT_MSG(bottom[6].value == 96, "\nscalar_cell, failed BOTTOM value");

        grid::partition_data<type> top = base_part.get_data(TOP).get();
        HPX_ASSERT_MSG(top.size() == x, "\nscalar_cell, failed TOP size");
        HPX_ASSERT_MSG(top.size_x() == x, "\nscalar_cell, failed TOP size_x");
        HPX_ASSERT_MSG(top.size_y() == 1, "\nscalar_cell, failed TOP size_y");
        HPX_ASSERT_MSG(top[8].value == 8, "\nscalar_cell, failed TOP value");

        grid::partition_data<type> bottom_left = base_part.get_data(BOTTOM_LEFT).get();
        HPX_ASSERT_MSG(bottom_left.size() == 1, "\nscalar_cell, failed BOTTOM_LEFT size");
        HPX_ASSERT_MSG(bottom_left.size_x() == 1, "\nscalar_cell, failed BOTTOM_LEFT size_x");
        HPX_ASSERT_MSG(bottom_left.size_y() == 1, "\nscalar_cell, failed BOTTOM_LEFT size_y");
        HPX_ASSERT_MSG(bottom_left[0].value == 99, "\nscalar_cell, failed BOTTOM_LEFT value");

        grid::partition_data<type> bottom_right = base_part.get_data(BOTTOM_RIGHT).get();
        HPX_ASSERT_MSG(bottom_right.size() == 1, "\nscalar_cell, failed BOTTOM_RIGHT size");
        HPX_ASSERT_MSG(bottom_right.size_x() == 1, "\nscalar_cell, failed BOTTOM_RIGHT size_x");
        HPX_ASSERT_MSG(bottom_right.size_y() == 1, "\nscalar_cell, failed BOTTOM_RIGHT size_y");
        HPX_ASSERT_MSG(bottom_right[0].value == 90, "\nscalar_cell, failed BOTTOM_RIGHT value");

        grid::partition_data<type> top_left = base_part.get_data(TOP_LEFT).get();
        HPX_ASSERT_MSG(top_left.size() == 1, "\nscalar_cell, failed TOP_LEFT size");
        HPX_ASSERT_MSG(top_left.size_x() == 1, "\nscalar_cell, failed TOP_LEFT size_x");
        HPX_ASSERT_MSG(top_left.size_y() == 1, "\nscalar_cell, failed TOP_LEFT size_y");
        HPX_ASSERT_MSG(top_left[0].value == 9, "\nscalar_cell, failed TOP_LEFT value");

        grid::partition_data<type> top_right = base_part.get_data(TOP_RIGHT).get();
        HPX_ASSERT_MSG(top_right.size() == 1, "\nscalar_cell, failed TOP_RIGHT size");
        HPX_ASSERT_MSG(top_right.size_x() == 1, "\nscalar_cell, failed TOP_RIGHT size_x");
        HPX_ASSERT_MSG(top_right.size_y() == 1, "\nscalar_cell, failed TOP_RIGHT size_y");
        HPX_ASSERT_MSG(top_right[0].value == 0, "\nscalar_cell, failed TOP_RIGHT value");
    }

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}
