#include "with_for_each.hpp"

#include <hpx/parallel/algorithm.hpp>
#include <hpx/parallel/algorithms/transform_reduce.hpp>

#include "util/helpers.hpp"
#include "stencils.hpp"

namespace computation {

void with_for_each::set_velocity_on_boundary(vector_data& uv_data, uint global_i, uint global_j,
                                                uint i_max, uint j_max)
{
    /*
    *@TODO: maybe create new vector_data here
    */
    uint size_x = uv_data.size_x();
    uint size_y = uv_data.size_y();

    bool is_left = (global_i == 0);
    bool is_right = (global_i + size_x > i_max);

    bool is_bottom = (global_j == 0);
    bool is_top = (global_j + size_y > j_max);

    uint start_i = (is_left ? 1 : 0);
    uint end_i = (is_right ? size_x - 1 : size_x);
    uint start_j = (is_bottom ? 1 : 0);
    uint end_j = (is_top ? size_y - 1 : size_y);

    auto range_i = boost::irange(start_i, end_i);
    auto range_j = boost::irange(start_j, end_j);

    vector_cell const topright_cell = uv_data.get_cell(size_x - 2, size_y - 2);

    std::vector<hpx::future<void> > futures;

    if (is_left)
    {
        futures.push_back(
            hpx::parallel::for_each(hpx::parallel::par(hpx::parallel::task), boost::begin(range_j), boost::end(range_j),
                [&uv_data](uint j)
                {
                    vector_cell& cell = uv_data.get_cell_ref(0, j);
                    vector_cell const cell2 = uv_data.get_cell(1, j);

                    cell.first = 0;
                    cell.second = -cell2.second;
                }
            )
        );
    }

    if (is_bottom)
    {
        futures.push_back(
            hpx::parallel::for_each(hpx::parallel::par(hpx::parallel::task), boost::begin(range_i), boost::end(range_i),
                [&uv_data](uint i)
                {
                    vector_cell& cell = uv_data.get_cell_ref(i, 0);
                    vector_cell const cell2 = uv_data.get_cell(i, 1);

                    cell.second = 0;
                    cell.first = -cell2.first;
                }
            )
        );
    }

    if (is_right)
    {
        futures.push_back(
            hpx::parallel::for_each(hpx::parallel::par(hpx::parallel::task), boost::begin(range_j), boost::end(range_j),
                [&uv_data, size_x](uint j)
                {
                    vector_cell& cell = uv_data.get_cell_ref(size_x - 2, j);
                    vector_cell& cell2 = uv_data.get_cell_ref(size_x - 1, j);

                    cell.first = 0;
                    cell2.second = -cell.second;
                }
            )
        );
    }

    if (is_top)
    {
        futures.push_back(
            hpx::parallel::for_each(hpx::parallel::par(hpx::parallel::task), boost::begin(range_i), boost::end(range_i),
                [&uv_data, size_y](uint i)
                {
                    vector_cell& cell = uv_data.get_cell_ref(i, size_y - 2);
                    vector_cell& cell2 = uv_data.get_cell_ref(i, size_y - 1);

                    cell.second = 0;
                    cell2.first = 2. - cell.first;
                }
            )
        );
    }

    hpx::wait_all(futures);

    if (is_right && is_top)
    {
        uv_data.get_cell_ref(size_x - 2, size_y - 1).first = 2. - topright_cell.first;
        uv_data.get_cell_ref(size_x - 1, size_y - 2).second = -topright_cell.second;
    }
}


}//computation

