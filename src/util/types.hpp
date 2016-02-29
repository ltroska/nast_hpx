#ifndef UTIL_TYPES_HPP
#define UTIL_TYPES_HPP

#define RealType double
#define uint std::size_t

#include <vector>
#include <tuple>


//forward declares
struct scalar_cell;
struct vector_cell;

namespace grid {
    template<typename T>
    struct partition_data;

    template<typename T>
    struct partition;
}

typedef grid::partition_data<scalar_cell> scalar_data;
typedef grid::partition_data<vector_cell> vector_data;
typedef grid::partition<scalar_cell> scalar_partition;
typedef grid::partition<vector_cell> vector_partition;
typedef std::vector<scalar_partition> scalar_grid_type;
typedef std::vector<vector_partition> vector_grid_type;
typedef std::vector<std::pair<uint, uint> > index_grid_type;

//directions for neighbor relations
enum direction
{
    LEFT = 0,
    TOP,
    BOTTOM_LEFT,
    BOTTOM_RIGHT,
    CENTER,
    TOP_LEFT,
    TOP_RIGHT,
    BOTTOM,
    RIGHT,
    NUM_DIRECTIONS
};

#endif
