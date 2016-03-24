#ifndef TEST_HELPERS_HPP
#define TEST_HELPERS_HPP

#include <random>

#include "grid/partition.hpp"
#include "computation/parameters.hpp"

class rand_double
{
public:
    rand_double(double low, double high)
    :r(std::bind(std::uniform_real_distribution<>(low,high),std::default_random_engine())){}

    double operator()(){ return r(); }

private:
    std::function<double()> r;
};

/*
* Makes the needed grids for testing the computation component
*/
struct grid_maker
{
    uint num_localities_x, num_localities_y;
    uint locality_id;

    computation::parameters p;

    grid_maker(uint num_localities_x_, uint num_localities_y_, uint locality_id_, computation::parameters params)
            : num_localities_x(num_localities_x_), num_localities_y(num_localities_y_), locality_id(locality_id_), p(params)
    {}

    void make_index_grid(index_grid_type& grid)
    {
        grid.resize(p.num_partitions_x * p.num_partitions_y);

        for (uint l = 0; l < p.num_partitions_y; l++)
            for (uint k = 0; k < p.num_partitions_x; k++)
            {
                grid[get_index(k, l)] =
                    std::pair<RealType, RealType>
                        (
                        (locality_id % num_localities_x) * (p.num_partitions_x - 2) * p.num_cells_per_partition_x + (k - 1) * p.num_cells_per_partition_x,
                        (locality_id / num_localities_x) * (p.num_partitions_y - 2) * p.num_cells_per_partition_y + (l - 1) * p.num_cells_per_partition_y
                        );
            }
    }

    void make_vector_grid(vector_grid_type& grid, RealType initial_value)
    {
        grid.resize(p.num_partitions_x * p.num_partitions_y);

        for (uint l = 0; l < p.num_partitions_y; l++)
            for (uint k = 0; k < p.num_partitions_x; k++)
            {
                grid[get_index(k, l)] = vector_partition(hpx::find_here(), p.num_cells_per_partition_x, p.num_cells_per_partition_y, initial_value);
            }
    }

    void make_scalar_grid(scalar_grid_type& grid, RealType initial_value)
    {
        grid.resize(p.num_partitions_x * p.num_partitions_y);

        for (uint l = 0; l < p.num_partitions_y; l++)
            for (uint k = 0; k < p.num_partitions_x; k++)
            {
                grid[get_index(k, l)] = scalar_partition(hpx::find_here(), p.num_cells_per_partition_x, p.num_cells_per_partition_y, initial_value);
            }
    }

    template<typename T>
    void make_random_grid(std::vector<grid::partition<T> >& grid, RealType lower_bound = -10, RealType upper_bound = 10)
    {
        rand_double rd{lower_bound, upper_bound};

        grid.resize(p.num_partitions_x * p.num_partitions_y);

        for (uint l = 0; l < p.num_partitions_y; l++)
            for (uint k = 0; k < p.num_partitions_x; k++)
            {
                grid::partition_data<T> data = grid::partition_data<T>(p.num_cells_per_partition_x, p.num_cells_per_partition_y);

                for (uint i = 0; i < p.num_cells_per_partition_x; i++)
                    for (uint j = 0; j < p.num_cells_per_partition_y; j++)
                    {
                        data.get_cell_ref(i, j) = T(rd());
                    }

                grid[get_index(k, l)] = grid::partition<T>(hpx::find_here(), data);
            }
    }

    template<typename T>
    void copy_grid(std::vector<grid::partition<T> > const & src_grid, std::vector<grid::partition<T> >& target_grid)
    {
        target_grid.resize(p.num_partitions_x * p.num_partitions_y);

        for (uint k = 0; k < p.num_partitions_x ; k++)
        {
            for (uint l = 0; l < p.num_partitions_y; l++)
            {
                grid::partition_data<T> base = grid::partition_data<T>(src_grid[get_index(k, l)].get_data(CENTER).get());

                target_grid[get_index(k, l)] = grid::partition<T>(hpx::find_here(), base);
            }
        }
    }

    void make_neighbor_test_grid(vector_grid_type& grid)
    {
        grid.resize(p.num_partitions_x * p.num_partitions_y);

        for (uint l = 0; l < p.num_partitions_y; l++)
            for (uint k = 0; k < p.num_partitions_x; k++)
            {
                vector_data data = vector_data(p.num_cells_per_partition_x, p.num_cells_per_partition_y);

                for (uint i = 0; i < p.num_cells_per_partition_x; i++)
                    for (uint j = 0; j < p.num_cells_per_partition_y; j++)
                    {
                        vector_cell& cell = data.get_cell_ref(i, j);

                        cell.first = k*p.num_cells_per_partition_x + i;
                        cell.second = l*p.num_cells_per_partition_y + j;

                    }

                grid[get_index(k, l)] = vector_partition(hpx::find_here(), data);
            }
    }

    private:
        uint get_index(uint k, uint l) { return l * p.num_partitions_x + k;}
};

template<typename T>
void print_grid(std::vector<grid::partition<T> > const& grid, computation::parameters p, uint start_k = 0, uint end_k = 0, uint start_l = 0, uint end_l = 0)
{
    std::vector<std::vector<grid::partition_data<T> > > data;

    data.resize(p.num_partitions_x);

    if (end_k == 0)
        end_k = p.num_partitions_x - 1;

    if (end_l == 0)
        end_l = p.num_partitions_y - 1;

    for (uint k = start_k; k <= end_k ; k++)
    {
        data[k].resize(p.num_partitions_y );
        for (uint l = start_l; l <= end_l; l++)
        {
            grid::partition_data<T> base = grid[l * p.num_partitions_x + k].get_data(CENTER).get();
            data[k][l] = grid::partition_data<T>(base);
        }
    }

    for (uint j = ((end_l > p.num_partitions_y -1) ? p.num_partitions_y - 1 : end_l); j <= p.num_partitions_y && j >= start_l; --j)
    {
        for (uint row = p.num_cells_per_partition_y - 1 ; row <= p.num_cells_per_partition_y; --row)
        {
            for (uint i = start_k; i <= end_k; ++i)
            {
                for (uint col = 0; col < p.num_cells_per_partition_x; col++)
                {
                    std::cout << data[i][j].get_cell(col, row) << " ";
                }
            }
            std::cout << std::endl;
        }
    }
}


std::string expected_string(RealType exp, RealType act)
{
    return "expected: " + std::to_string(exp) + " got: " + std::to_string(act);
}

std::string expected_string(vector_cell cell, vector_cell neighbor)
{
    return "cell: i=" + std::to_string(cell.first) + " j=" + std::to_string(cell.second)
            + "| neighbor: i=" + std::to_string(neighbor.first) + " j=" + std::to_string(neighbor.second);
}

#endif
