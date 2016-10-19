#ifndef NAST_HPX_GRID_PARTITION_DATA_HPP
#define NAST_HPX_GRID_PARTITION_DATA_HPP

#include <hpx/runtime/serialization/serialize.hpp>

#include "util/defines.hpp"

namespace nast_hpx { namespace grid {

/// This enum represents all directions needed for 2D neighbor relations
enum direction
{
    LEFT = 0,
    TOP,
    BOTTOM_RIGHT,
    TOP_LEFT,
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

    partition_data(std::size_t size_x, std::size_t size_y)
    : data_(size_x * size_y, T()),
      size_x_(size_x),
      size_y_(size_y),
      size_(size_x * size_y)
    {}

    void resize(std::size_t size_x, std::size_t size_y, T val = T())
    {
        size_x_ = size_x;
        size_y_ = size_y;
        size_ = size_x * size_y;
        data_.resize(size_x_ * size_y_, val);
    }

    void clear(T val = T())
    {
        for (auto& elem : data_)
            elem = val;
    }

    inline T operator[](std::size_t idx) const { return data_[idx];}
    inline T& operator[](std::size_t idx) { return data_[idx];}

    inline T& operator()(unsigned idx, unsigned idy)
    {return data_[idy * size_x_ + idx];}

    inline T const& operator()(unsigned idx, unsigned idy) const
    {return data_[idy * size_x_ + idx];}

    typename std::vector<T>::iterator begin() { return data_.begin(); }
    typename std::vector<T>::iterator end() { return data_.end(); }

    friend std::ostream& operator<<(std::ostream& os,
        partition_data<T> const& data)
    {
        os << data.size_x_ << " " << data.size_y_ << " " << data.size_ << "\n";
        for (std::size_t j = data.size_y_ - 1; j < data.size_y_; --j)
        {
            for (std::size_t i = 0; i < data.size_x_; ++i)
                os << data(i, j) << " ";

            os << "\n";
        }

        return os;
    }

    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned version)
    {
        ar & size_x_ & size_y_ & size_ & data_;
    }

    std::vector<T> data_;
    std::size_t size_x_;
    std::size_t size_y_;
    std::size_t size_;
};

}//namespace grid
}

#endif
