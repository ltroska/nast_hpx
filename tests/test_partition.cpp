#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include "grid/partition.hpp"

int hpx_main(int argc, char* argv[])
{
    std::size_t x = 5;
    std::size_t y = 5;
    std::size_t size = x*y;

    //test empty
    grid::partition empty_part(hpx::find_here(), x, y);
    grid::partition_data<> empty_pdata = empty_part.get_data(CENTER).get();

    for (int i = 0; i < size; i++)
        HPX_ASSERT(empty_pdata[i] == 0);

    //test initial value
    grid::partition initial_part(hpx::find_here(), x, y, 4);
    grid::partition_data<> initial_pdata = initial_part.get_data(CENTER).get();

    for (int i = 0; i < size; i++)
        HPX_ASSERT(initial_pdata[i] == 4);

    //test base data
    grid::partition_data<> base_pdata(5, 5);
    for (int i = 0; i < size; i++)
            base_pdata[i] = i;

    grid::partition base_part(hpx::find_here(), base_pdata);
    grid::partition_data<> base_recv_pdata = base_part.get_data(CENTER).get();
    for (int i = 0; i < size; i++)
            HPX_ASSERT(base_recv_pdata[i] == base_pdata[i]);

    //test directional
    grid::partition_data<> left = base_part.get_data(LEFT).get();
    HPX_ASSERT(left.size() == 5);
    HPX_ASSERT(left.size_x() == 1);
    HPX_ASSERT(left.size_y() == 5);
    HPX_ASSERT(left[1] == 9);

    grid::partition_data<> right = base_part.get_data(RIGHT).get();
    HPX_ASSERT(right.size() == 5);
    HPX_ASSERT(right.size_x() == 1);
    HPX_ASSERT(right.size_y() == 5);
    HPX_ASSERT(right[2] == 10);

    grid::partition_data<> bottom = base_part.get_data(BOTTOM).get();
    HPX_ASSERT(bottom.size() == 5);
    HPX_ASSERT(bottom.size_x() == 5);
    HPX_ASSERT(bottom.size_y() == 1);
    HPX_ASSERT(bottom[3] == 23);

    grid::partition_data<> top = base_part.get_data(TOP).get();
    HPX_ASSERT(top.size() == 5);
    HPX_ASSERT(top.size_x() == 5);
    HPX_ASSERT(top.size_y() == 1);
    HPX_ASSERT(top[3] == 3);

    grid::partition_data<> bottom_left = base_part.get_data(BOTTOM_LEFT).get();
    HPX_ASSERT(bottom_left.size() == 1);
    HPX_ASSERT(bottom_left.size_x() == 1);
    HPX_ASSERT(bottom_left.size_y() == 1);
    HPX_ASSERT(bottom_left[0] == 24);

    grid::partition_data<> bottom_right = base_part.get_data(BOTTOM_RIGHT).get();
    HPX_ASSERT(bottom_right.size() == 1);
    HPX_ASSERT(bottom_right.size_x() == 1);
    HPX_ASSERT(bottom_right.size_y() == 1);
    HPX_ASSERT(bottom_right[0] == 20);

    grid::partition_data<> top_left = base_part.get_data(TOP_LEFT).get();
    HPX_ASSERT(top_left.size() == 1);
    HPX_ASSERT(top_left.size_x() == 1);
    HPX_ASSERT(top_left.size_y() == 1);
    HPX_ASSERT(top_left[0] == 4);

    grid::partition_data<> top_right = base_part.get_data(TOP_RIGHT).get();
    HPX_ASSERT(top_right.size() == 1);
    HPX_ASSERT(top_right.size_x() == 1);
    HPX_ASSERT(top_right.size_y() == 1);
    HPX_ASSERT(top_right[0] == 0);

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}
