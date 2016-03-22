#ifndef GRID_TYPES_HPP
#define GRID_TYPES_HPP

#include <vector>

#include "cell.hpp"
#include "partition_data.hpp"
#include "partition.hpp"

typedef grid::scalar_cell scalar_cell;
typedef grid::vector_cell vector_cell;
typedef grid::partition_data<scalar_cell> scalar_data;
typedef grid::partition_data<vector_cell> vector_data;
typedef grid::partition<scalar_cell> scalar_partition;
typedef grid::partition<vector_cell> vector_partition;
typedef std::vector<scalar_partition> scalar_grid_type;
typedef std::vector<vector_partition> vector_grid_type;
typedef std::vector<std::pair<uint, uint> > index_grid_type;

#endif
