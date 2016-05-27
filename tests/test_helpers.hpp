#ifndef TEST_HELPERS_HPP
#define TEST_HELPERS_HPP

#include <random>

#include "grid/types.hpp"
#include "util/helpers.hpp"

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
    
    uint cells_x, cells_y, partitions_x, partitions_y;

    grid_maker(uint num_localities_x_, uint num_localities_y_,
                uint locality_id_, uint cells_x_, uint cells_y_,
                uint partitions_x_, uint partitions_y_)
            : num_localities_x(num_localities_x_),
                num_localities_y(num_localities_y_),
                locality_id(locality_id_),
                cells_x(cells_x_),
                cells_y(cells_y_),
                partitions_x(partitions_x_),
                partitions_y(partitions_y_)
    {}

    void make_index_grid(index_grid_type& grid)
    {
        grid.resize(partitions_x * partitions_y);

        for (uint l = 0; l < partitions_y; l++)
            for (uint k = 0; k < partitions_x; k++)
            {
                grid[get_index(k, l)] =
                    std::pair<RealType, RealType>
                        (
                        (locality_id % num_localities_x) * (partitions_x - 2) * cells_x + (k - 1) * cells_x,
                        (locality_id / num_localities_x) * (partitions_y - 2) * cells_y + (l - 1) * cells_y
                        );
            }
    }
    
    void make_flag_grid(flag_grid_type& flag_grid, uint i_max, uint j_max)
    {
        flag_grid.resize(partitions_x * partitions_y);
        flag_grid[0].resize(cells_x * cells_y);
        
        for (uint l = 0; l < partitions_y; ++l)
            for (uint k = 0; k < partitions_x; ++k)
                for (uint j = 0; j < cells_y; ++j)
                    for (uint i = 0; i < cells_x; i++)
                    {
                        uint global_i = (locality_id % num_localities_x) * (partitions_x - 2) * cells_x + (k - 1) * cells_x + i;
                        uint global_j = (locality_id / num_localities_x) * (partitions_y - 2) * cells_y + (l - 1) * cells_y + j;
                        
                        if (in_range(0, 0, 1, j_max, global_i, global_j))
                            flag_grid[l * partitions_x + k].push_back(std::bitset<5>("00111"));
                        else if (in_range(i_max + 1, i_max + 1, 1, j_max, global_i, global_j))
                            flag_grid[l * partitions_x + k].push_back(std::bitset<5>("01011"));
                        else if (in_range(1, i_max, 0, 0, global_i, global_j))
                            flag_grid[l * partitions_x + k].push_back(std::bitset<5>("01110"));
                        else if (in_range(1, i_max, j_max + 1, j_max + 1, global_i, global_j))
                            flag_grid[l * partitions_x + k].push_back(std::bitset<5>("01101"));
                        else if (in_range(1, 1, 2, j_max - 1, global_i, global_j))
                            flag_grid[l * partitions_x + k].push_back(std::bitset<5>("11011"));
                        else if (in_range(i_max, i_max, 2, j_max - 1, global_i, global_j))
                            flag_grid[l * partitions_x + k].push_back(std::bitset<5>("10111"));
                        else if (in_range(2, i_max - 1, 1, 1, global_i, global_j))
                            flag_grid[l * partitions_x + k].push_back(std::bitset<5>("11101"));
                        else if (in_range(2, i_max - 1, j_max, j_max, global_i, global_j))
                            flag_grid[l * partitions_x + k].push_back(std::bitset<5>("11110"));
                        else if (in_range(1, 1, 1, 1, global_i, global_j))
                            flag_grid[l * partitions_x + k].push_back(std::bitset<5>("11001"));
                        else if (in_range(1, 1, j_max, j_max, global_i, global_j))
                            flag_grid[l * partitions_x + k].push_back(std::bitset<5>("11010"));
                        else if (in_range(i_max, i_max, 1, 1, global_i, global_j))
                            flag_grid[l * partitions_x + k].push_back(std::bitset<5>("10101"));
                        else if (in_range(i_max, i_max, j_max, j_max, global_i, global_j))
                            flag_grid[l * partitions_x + k].push_back(std::bitset<5>("10110"));
                        else
                            flag_grid[l * partitions_x + k].push_back(std::bitset<5>("11111"));    
                    }
    }

    void make_vector_grid(vector_grid_type& grid, RealType initial_value)
    {
        grid.resize(partitions_x * partitions_y);

        for (uint l = 0; l < partitions_y; l++)
            for (uint k = 0; k < partitions_x; k++)
            {
                grid[get_index(k, l)] = vector_partition(hpx::find_here(), cells_x, cells_y, initial_value);
            }
    }

    void make_scalar_grid(scalar_grid_type& grid, RealType initial_value)
    {
        grid.resize(partitions_x * partitions_y);

        for (uint l = 0; l < partitions_y; l++)
            for (uint k = 0; k < partitions_x; k++)
            {
                grid[get_index(k, l)] = scalar_partition(hpx::find_here(), cells_x, cells_y, initial_value);
            }
    }

    template<typename T>
    void make_random_grid(std::vector<grid::partition<T> >& grid, RealType lower_bound = -10, RealType upper_bound = 10)
    {
        rand_double rd{lower_bound, upper_bound};

        grid.resize(partitions_x * partitions_y);

        for (uint l = 0; l < partitions_y; l++)
            for (uint k = 0; k < partitions_x; k++)
            {
                grid::partition_data<T> data = grid::partition_data<T>(cells_x, cells_y);

                for (uint i = 0; i < cells_x; i++)
                    for (uint j = 0; j < cells_y; j++)
                    {
                        data(i, j) = T(rd());
                    }

                grid[get_index(k, l)] = grid::partition<T>(hpx::find_here(), data);
            }
    }

    template<typename T>
    void copy_grid(std::vector<grid::partition<T> > const & src_grid, std::vector<grid::partition<T> >& target_grid)
    {
        target_grid.resize(partitions_x * partitions_y);

        for (uint k = 0; k < partitions_x ; k++)
        {
            for (uint l = 0; l < partitions_y; l++)
            {
                auto base_data = src_grid[get_index(k, l)].get_data(CENTER).get();
                
                grid::partition_data<T> base(base_data.size_x(), base_data.size_y());
                
                for (uint idx = 0; idx < base_data.size(); ++idx)
                    base[idx] = base_data[idx];
                
                target_grid[get_index(k, l)] = grid::partition<T>(hpx::find_here(), base);  
            }
        }
    }

    void make_neighbor_test_grid(vector_grid_type& grid)
    {
        grid.resize(partitions_x * partitions_y);

        for (uint l = 0; l < partitions_y; l++)
            for (uint k = 0; k < partitions_x; k++)
            {
                vector_data data = vector_data(cells_x, cells_y);

                for (uint i = 0; i < cells_x; i++)
                    for (uint j = 0; j < cells_y; j++)
                    {
                        vector_cell& cell = data(i, j);

                        cell.first = k*cells_x + i;
                        cell.second = l*cells_y + j;

                    }

                grid[get_index(k, l)] = vector_partition(hpx::find_here(), data);
            }
    }

    private:
        uint get_index(uint k, uint l) { return l * partitions_x + k;}
};

template<typename T>
void print_grid(std::vector<grid::partition<T> > const& grid, uint partitions_x,
    uint start_k = 0, uint end_k = 0, uint start_l = 0, uint end_l = 0)
{
    std::vector<std::vector<grid::partition_data<T> > > data;
    
    uint partitions_y = grid.size()/partitions_x;
    
    uint cells_x, cells_y;

    data.resize(partitions_x);

    if (end_k == 0)
        end_k = partitions_x - 1;

    if (end_l == 0)
        end_l = partitions_y - 1;

    for (uint k = start_k; k <= end_k ; k++)
    {
        data[k].resize(partitions_y );
        for (uint l = start_l; l <= end_l; l++)
        {
            grid::partition_data<T> base = grid[l * partitions_x + k].get_data(CENTER).get();
            data[k][l] = grid::partition_data<T>(base);
            
            if (k == 1 && l == 1)
            {
                cells_x = base.size_x();
                cells_y = base.size_y();
            }
        }
    }

    for (uint j = ((end_l > partitions_y -1) ? partitions_y - 1 : end_l); j <= partitions_y && j >= start_l; --j)
    {
        for (uint row = cells_y - 1 ; row <= cells_y; --row)
        {
            for (uint i = start_k; i <= end_k; ++i)
            {
                for (uint col = 0; col < cells_x; col++)
                {
                    std::cout << data[i][j](col, row) << " ";
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
