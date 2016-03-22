#ifndef UTIL_HELPERS_HPP
#define UTIL_HELPERS_HPP

#include "typedefs.hpp"

//finding id of neighboring localities
inline const uint get_neighbor_id(uint id, direction dir, uint num_localities)
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
            return ((id-1)/res_x == id/res_x && id-1 < num_localities) ? id-1 : num_localities;

        case RIGHT:
            return ((id+1)/res_x == id/res_x && id+1 < num_localities) ? id+1 : num_localities;

        case TOP:
            return ((id+res_x) < num_localities && (id+res_x)/res_x == id/res_x +1) ? id+res_x : num_localities;

        case BOTTOM:
            return ((id-res_x) < num_localities && (id-res_x)/res_x == id/res_x -1) ? id-res_x : num_localities;

        case TOP_LEFT:
            return ((id+res_x-1) < num_localities && (id+res_x-1)/res_x == id/res_x +1) ? id+res_x-1 : num_localities;

        case TOP_RIGHT:
            return ((id+res_x+1) < num_localities && (id+res_x+1)/res_x == id/res_x+1) ? id+res_x+1 : num_localities;

        case BOTTOM_LEFT:
            return ((id-res_x-1) < num_localities && (id-res_x-1)/res_x == id/res_x-1) ? id-res_x-1 : num_localities;

        case BOTTOM_RIGHT:
            return ((id-res_x+1) < num_localities && (id-res_x+1)/res_x == id/res_x-1) ? id-res_x+1 : num_localities;
        case CENTER:
            return id;
        default:
            return num_localities;
    }
}

//check if point is in global id range
inline bool in_range(uint start_i, uint end_i, uint start_j, uint end_j, uint i,  uint j)
{
    return !(i < start_i || i > end_i || j < start_j || j > end_j);
}

template<typename T>
const T get_left_neighbor(grid::partition_data<T> const& center, grid::partition_data<T> const& left, uint i, uint j)
{
    if (i > 0)
        return center.get_cell(i-1, j);
    else
        return left.get_cell(j, 0);
}

template<typename T>
const T get_right_neighbor(grid::partition_data<T> const& center, grid::partition_data<T> const& right, uint i, uint j)
{
    if (i+1 < center.size_x())
        return center.get_cell(i+1, j);
    else
        return right.get_cell(j, 0);
}

template<typename T>
const T get_bottom_neighbor(grid::partition_data<T> const& center, grid::partition_data<T> const& bottom, uint i, uint j)
{
    if (j > 0)
        return center.get_cell(i, j-1);
    else
        return bottom.get_cell(i, 0);
}

template<typename T>
const T get_top_neighbor(grid::partition_data<T> const& center, grid::partition_data<T> const& top, uint i, uint j)
{
    if (j+1 < center.size_y())
        return center.get_cell(i, j+1);
    else
        return top.get_cell(i, 0);
}

template<typename T>
const T get_neighbor_cell(grid::partition_data<T> const& center, grid::partition_data<T> const& left, grid::partition_data<T> const& right,
                           grid::partition_data<T> const& bottom, grid::partition_data<T> const& top, grid::partition_data<T> const& bottomleft,
                           grid::partition_data<T> const& bottomright, grid::partition_data<T> const& topleft, grid::partition_data<T> const& topright,
                           uint i, uint j, direction dir)
{

    uint size_x = center.size_x();
    uint size_y = center.size_y();

    switch (dir)
    {
        case LEFT:
            if (i > 0)
                return center.get_cell(i-1, j);
            else
                return left.get_cell(j, 0);

        case RIGHT:
            if (i+1 < size_x)
                return center.get_cell(i+1, j);
            else
                return right.get_cell(j, 0);

        case BOTTOM:
            if (j > 0)
                return center.get_cell(i, j-1);
            else
                return bottom.get_cell(i, 0);

        case TOP:
            if (j+1 < size_y)
                return center.get_cell(i, j+1);
            else
                return top.get_cell(i, 0);

        case BOTTOM_LEFT:
            if (i > 0 && j > 0)
                return center.get_cell(i-1, j-1);
            else if (i > 0)
                return bottom.get_cell(i-1, 0);
            else if (j > 0)
                return left.get_cell(j-1, 0);
            else
                return bottomleft.get_cell(0, 0);

        case BOTTOM_RIGHT:
            if (i+1 < size_x && j > 0)
                return center.get_cell(i+1, j-1);
            else if (i+1 < size_x)
                return bottom.get_cell(i+1, 0);
            else if (j > 0)
                return right.get_cell(j-1, 0);
            else
                return bottomright.get_cell(0, 0);

        case TOP_LEFT:
            if (i > 0 && j+1 < size_y)
                return center.get_cell(i-1, j+1);
            else if (i > 0)
                return top.get_cell(i-1, 0);
            else if (j+1 < size_y)
                return left.get_cell(j+1, 0);
            else
                return topleft.get_cell(0, 0);

        case TOP_RIGHT:
            if (i+1 < size_x && j+1 < size_y)
                return center.get_cell(i+1, j+1);
            else if (i+1 < size_x)
                return top.get_cell(i+1, 0);
            else if (j+1 < size_y)
                return right.get_cell(j+1, 0);
            else
                return topright.get_cell(0, 0);

        default:
            return center.get_cell(i, j);
    }
}
#endif
