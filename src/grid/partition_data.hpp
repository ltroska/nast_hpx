#ifndef NAST_HPX_GRID_PARTITION_DATA_HPP
#define NAST_HPX_GRID_PARTITION_DATA_HPP

#include <hpx/runtime/serialization/serialize.hpp>

#include "util/typedefs.hpp"

namespace nast_hpx { namespace grid {

/// This enum represents all directions needed for 2D neighbor relations
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
    
/// This class represents a block of a grid.
template<typename T = Real>
struct partition_data
{
public:

    partition_data()
    : size_x_(0),
      size_y_(0),
      size_(0)
    {}

    partition_data(uint size_x, uint size_y)
    : data_(size_x * size_y, 0),
      size_x_(size_x),
      size_y_(size_y),
      size_(size_x * size_y),
      act_size_x_(size_x_ - 2),
      act_size_y_(size_y_ - 2),
      act_size_(act_size_x_ * act_size_y_)
    {}

    partition_data(uint size_x, uint size_y, T initial_value)
    : data_(size_x * size_y, initial_value),
      size_x_(size_x),
      size_y_(size_y),
      size_(size_x * size_y),
      act_size_x_(size_x_ - 2),
      act_size_y_(size_y_ - 2),
      act_size_(act_size_x_ * act_size_y_)
    {}

    inline T operator[](uint idx) const { return data_[idx];}
    inline T& operator[](uint idx) { return data_[idx];}
    
    inline T& operator()(unsigned idx, unsigned idy)
    {return data_[idy * size_x_ + idx];}
    
    inline T const& operator()(unsigned idx, unsigned idy) const
    {return data_[idy * size_x_ + idx];}

    T* begin() { return data_.begin(); }
    T* end() { return data_.end(); }

    friend std::ostream& operator<<(std::ostream& os,
        partition_data<T> const& data)
    {
        for (uint j = data.size_y_ - 1; j < data.size_y_; --j)
        {
            for (uint i = 0; i < data.size_x_; ++i)
                os << data(i, j) << " ";

            os << "\n";
        }

        return os;
    }
    
    std::vector<T> data_;
    uint size_x_;
    uint size_y_;
    uint size_;
    
    uint act_size_x_;
    uint act_size_y_;
    uint act_size_;
};

}//namespace grid
}

#endif
