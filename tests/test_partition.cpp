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
            HPX_ASSERT(empty_pdata[i] == 0);

        //test initial value
        grid::partition<type> initial_part(hpx::find_here(), x, y, 4);
        grid::partition_data<type> initial_pdata = initial_part.get_data(CENTER).get();

        for (int i = 0; i < size; i++)
            HPX_ASSERT(initial_pdata[i] == 4);

        //test base data
        grid::partition_data<type> base_pdata(5, 5);
        for (int i = 0; i < size; i++)
                base_pdata[i] = i;

        grid::partition<type> base_part(hpx::find_here(), base_pdata);
        grid::partition_data<type> base_recv_pdata = base_part.get_data(CENTER).get();
        for (int i = 0; i < size; i++)
                HPX_ASSERT(base_recv_pdata[i] == base_pdata[i]);

        //test directional
        grid::partition_data<type> left = base_part.get_data(LEFT).get();
        HPX_ASSERT(left.size() == y);
        HPX_ASSERT(left.size_x() == 1);
        HPX_ASSERT(left.size_y() == y);
        HPX_ASSERT(left[1] == 9);

        grid::partition_data<type> right = base_part.get_data(RIGHT).get();
        HPX_ASSERT(right.size() == y);
        HPX_ASSERT(right.size_x() == 1);
        HPX_ASSERT(right.size_y() == y);
        HPX_ASSERT(right[2] == 10);

        grid::partition_data<type> bottom = base_part.get_data(BOTTOM).get();
        HPX_ASSERT(bottom.size() == x);
        HPX_ASSERT(bottom.size_x() == x);
        HPX_ASSERT(bottom.size_y() == 1);
        HPX_ASSERT(bottom[3] == 23);

        grid::partition_data<type> top = base_part.get_data(TOP).get();
        HPX_ASSERT(top.size() == x);
        HPX_ASSERT(top.size_x() == x);
        HPX_ASSERT(top.size_y() == 1);
        HPX_ASSERT(top[3] == 3);

        grid::partition_data<type> bottom_left = base_part.get_data(BOTTOM_LEFT).get();
        HPX_ASSERT(bottom_left.size() == 1);
        HPX_ASSERT(bottom_left.size_x() == 1);
        HPX_ASSERT(bottom_left.size_y() == 1);
        HPX_ASSERT(bottom_left[0] == 24);

        grid::partition_data<type> bottom_right = base_part.get_data(BOTTOM_RIGHT).get();
        HPX_ASSERT(bottom_right.size() == 1);
        HPX_ASSERT(bottom_right.size_x() == 1);
        HPX_ASSERT(bottom_right.size_y() == 1);
        HPX_ASSERT(bottom_right[0] == 20);

        grid::partition_data<type> top_left = base_part.get_data(TOP_LEFT).get();
        HPX_ASSERT(top_left.size() == 1);
        HPX_ASSERT(top_left.size_x() == 1);
        HPX_ASSERT(top_left.size_y() == 1);
        HPX_ASSERT(top_left[0] == 4);

        grid::partition_data<type> top_right = base_part.get_data(TOP_RIGHT).get();
        HPX_ASSERT(top_right.size() == 1);
        HPX_ASSERT(top_right.size_x() == 1);
        HPX_ASSERT(top_right.size_y() == 1);
        HPX_ASSERT(top_right[0] == 0);
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
            HPX_ASSERT(empty_pdata[i].value == 0);

        //test initial value
        grid::partition<type> initial_part(hpx::find_here(), x, y, 4);
        grid::partition_data<type> initial_pdata = initial_part.get_data(CENTER).get();

        for (int i = 0; i < size; i++)
            HPX_ASSERT(initial_pdata[i].value == 4);

        //test base data
        grid::partition_data<type> base_pdata(x, y);
        for (int i = 0; i < size; i++)
                base_pdata[i].value = i;

        grid::partition<type> base_part(hpx::find_here(), base_pdata);
        grid::partition_data<type> base_recv_pdata = base_part.get_data(CENTER).get();
        for (int i = 0; i < size; i++)
                HPX_ASSERT(base_recv_pdata[i].value == base_pdata[i].value);

        //test directional
        grid::partition_data<type> left = base_part.get_data(LEFT).get();
        HPX_ASSERT(left.size() == y);
        HPX_ASSERT(left.size_x() == 1);
        HPX_ASSERT(left.size_y() == y);
        HPX_ASSERT(left[3].value == 39);

        grid::partition_data<type> right = base_part.get_data(RIGHT).get();
        HPX_ASSERT(right.size() == y);
        HPX_ASSERT(right.size_x() == 1);
        HPX_ASSERT(right.size_y() == y);
        HPX_ASSERT(right[4].value == 40);

        grid::partition_data<type> bottom = base_part.get_data(BOTTOM).get();
        HPX_ASSERT(bottom.size() == x);
        HPX_ASSERT(bottom.size_x() == x);
        HPX_ASSERT(bottom.size_y() == 1);
        HPX_ASSERT(bottom[6].value == 96);

        grid::partition_data<type> top = base_part.get_data(TOP).get();
        HPX_ASSERT(top.size() == x);
        HPX_ASSERT(top.size_x() == x);
        HPX_ASSERT(top.size_y() == 1);
        HPX_ASSERT(top[8].value == 8);

        grid::partition_data<type> bottom_left = base_part.get_data(BOTTOM_LEFT).get();
        HPX_ASSERT(bottom_left.size() == 1);
        HPX_ASSERT(bottom_left.size_x() == 1);
        HPX_ASSERT(bottom_left.size_y() == 1);
        HPX_ASSERT(bottom_left[0].value == 99);

        grid::partition_data<type> bottom_right = base_part.get_data(BOTTOM_RIGHT).get();
        HPX_ASSERT(bottom_right.size() == 1);
        HPX_ASSERT(bottom_right.size_x() == 1);
        HPX_ASSERT(bottom_right.size_y() == 1);
        HPX_ASSERT(bottom_right[0].value == 90);

        grid::partition_data<type> top_left = base_part.get_data(TOP_LEFT).get();
        HPX_ASSERT(top_left.size() == 1);
        HPX_ASSERT(top_left.size_x() == 1);
        HPX_ASSERT(top_left.size_y() == 1);
        HPX_ASSERT(top_left[0].value == 9);

        grid::partition_data<type> top_right = base_part.get_data(TOP_RIGHT).get();
        HPX_ASSERT(top_right.size() == 1);
        HPX_ASSERT(top_right.size_x() == 1);
        HPX_ASSERT(top_right.size_y() == 1);
        HPX_ASSERT(top_right[0].value == 0);
    }

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}
