#ifndef UTIL_HELPERS_HPP
#define UTIL_HELPERS_HPP

#include "typedefs.hpp"

/// Method that finds the id of the neighboring locality
uint inline get_neighbor_id(uint id, direction dir, uint num_localities)
{
    uint res_x, res_y;
    if (num_localities == 2)
    {
        res_x = 2;
        res_y = 1;
    }
    else
    {
        res_x = static_cast<uint>(sqrt(num_localities));
        res_y = res_x;
    }

    switch (dir)
    {
        case LEFT:
            return ((id-1)/res_x == id/res_x && id-1 < num_localities)
                ? id-1 : num_localities;

        case RIGHT:
            return ((id+1)/res_x == id/res_x && id+1 < num_localities)
                ? id+1 : num_localities;

        case TOP:
            return ((id+res_x) < num_localities
                        && (id+res_x)/res_x == id/res_x +1)
                ? id+res_x : num_localities;

        case BOTTOM:
            return ((id-res_x) < num_localities
                        && (id-res_x)/res_x == id/res_x -1)
                ? id-res_x : num_localities;

        case TOP_LEFT:
            return ((id+res_x-1) < num_localities
                        && (id+res_x-1)/res_x == id/res_x +1)
                ? id+res_x-1 : num_localities;

        case TOP_RIGHT:
            return ((id+res_x+1) < num_localities
                        && (id+res_x+1)/res_x == id/res_x+1)
                ? id+res_x+1 : num_localities;

        case BOTTOM_LEFT:
            return ((id-res_x-1) < num_localities
                        && (id-res_x-1)/res_x == id/res_x-1)
                ? id-res_x-1 : num_localities;

        case BOTTOM_RIGHT:
            return ((id-res_x+1) < num_localities
                        && (id-res_x+1)/res_x == id/res_x-1)
                ? id-res_x+1 : num_localities;
            
        case CENTER:
            return id;
        default:
            return num_localities;
    }
}

/// Method to check if cell is in global index range
bool inline in_range(uint start_i, uint end_i, uint start_j, uint end_j, uint i,
                        uint j)
{
    return !(i < start_i || i > end_i || j < start_j || j > end_j);
}

/// Method extracts the left neighbor of the given cell from the provided
/// partitions.
template<typename T>
inline T get_left_neighbor(grid::partition_data<T> const& center,
                        grid::partition_data<T> const& left, uint i, uint j)
{
    if (i > 0)
        return center(i-1, j);
    else
        return left[j];
}

/// Method extracts the right neighbor of the given cell from the provided
/// partitions
template<typename T>
inline T get_right_neighbor(grid::partition_data<T> const& center,
                        grid::partition_data<T> const& right, uint i, uint j)
{
    if (i+1 < center.size_x())
        return center(i+1, j);
    else
        return right[j];
}

/// Method extracts the bottom neighbor of the given cell from the provided
/// partitions
template<typename T>
inline T get_bottom_neighbor(grid::partition_data<T> const& center,
                        grid::partition_data<T> const& bottom, uint i, uint j)
{
    if (j > 0)
        return center(i, j-1);
    else
        return bottom[i];
}

/// Method extracts the top neighbor of the given cell from the provided
/// partitions
template<typename T>
inline T get_top_neighbor(grid::partition_data<T> const& center,
                    grid::partition_data<T> const& top, uint i, uint j)
{
    if (j+1 < center.size_y())
        return center(i, j+1);
    else
        return top[i];
}

/// Method extracts the bottom right neighbor of the given cell from the
/// provided partitions
template<typename T>
T get_bottomright_neighbor(grid::partition_data<T> const& center, 
        grid::partition_data<T> const& bottom,
        grid::partition_data<T> const& right,
        grid::partition_data<T> const& bottomright, uint i, uint j)
{
    if (i+1 < center.size_x() && j > 0)
        return center(i+1, j-1);
    else if (i+1 < center.size_x())
        return bottom(i+1, 0);
    else if (j > 0)
        return right(j-1, 0);
    else
        return bottomright(0, 0);
}

/// Method extracts the top left neighbor of the given cell from the
/// provided partitions
template<typename T>
T get_topleft_neighbor(grid::partition_data<T> const& center,
        grid::partition_data<T> const& top,
        grid::partition_data<T> const& left,
        grid::partition_data<T> const& topleft, uint i, uint j)
{
    if (i > 0 && j+1 < center.size_y())
        return center(i-1, j+1);
    else if (i > 0)
        return top(i-1, 0);
    else if (j+1 < center.size_y())
        return left(j+1, 0);
    else
        return topleft(0, 0);
}
#endif
