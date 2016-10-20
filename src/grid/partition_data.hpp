#ifndef NAST_HPX_GRID_PARTITION_DATA_HPP_
#define NAST_HPX_GRID_PARTITION_DATA_HPP_

#include "util/defines.hpp"

#include "util/hpx_wrap.hpp"

namespace nast_hpx { namespace grid {

/// This class represents a block of a grid.
template<typename T = double>
struct partition_data
{
public:

    partition_data()
    : size_x_(0),
      size_y_(0),
      size_z_(0),
      size_(0)
    {}

    partition_data(std::size_t size_x, std::size_t size_y, std::size_t size_z, T val = T())
    : data_(size_x * size_y * size_z, val),
      size_x_(size_x),
      size_y_(size_y),
      size_z_(size_z),
      size_(size_x * size_y * size_z)
    {}

    void resize(std::size_t size_x, std::size_t size_y, std::size_t size_z, T val = T())
    {
        size_x_ = size_x;
        size_y_ = size_y;
        size_z_ = size_z;
        size_ = size_x * size_y * size_z;
        data_.resize(size_x_ * size_y_ * size_z_, val);
    }

    void clear(T val = T())
    {
        for (auto& elem : data_)
            elem = val;
    }

    inline T operator[](std::size_t idx) const { return data_[idx];}
    inline T& operator[](std::size_t idx) { return data_[idx];}

    inline T& operator()(std::size_t idx, std::size_t idy, std::size_t idz)
    {return data_[idz * size_y_ * size_x_ + idy * size_x_ + idx];}

    inline T const& operator()(std::size_t idx, std::size_t idy, std::size_t idz) const
    {return data_[idz * size_y_ * size_x_ + idy * size_x_ + idx];}

    typename std::vector<T>::iterator begin() { return data_.begin(); }
    typename std::vector<T>::iterator end() { return data_.end(); }

    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned version)
    {
        ar & size_x_ & size_y_ & size_z_ & size_ & data_;
    }

    std::vector<T> data_;
    std::size_t size_x_;
    std::size_t size_y_;
    std::size_t size_z_;
    std::size_t size_;
};

}//namespace grid
}

#endif
