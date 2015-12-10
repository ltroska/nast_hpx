#ifndef GRID_PARTITION_DATA_HPP
#define GRID_PARTITION_DATA_HPP

#include <hpx/runtime/serialization/serialize.hpp>
#include <boost/shared_array.hpp>

#include "internal/types.hpp"
#include "internal/utils.hpp"
#include "cell.hpp"

namespace grid {

//flag for which data we are requesting
enum partition_type
{
    left_partition, center_partition, right_partition, top_partition, bottom_partition, top_left_partition, top_right_partition,
    bottom_left_partition, bottom_right_partition
};

//partition_data
template<typename T>
struct partition_data
{
private:
    typedef hpx::serialization::serialize_buffer<T> buffer_type;

public:
    partition_data()
    : size_x_(0),
      size_y_(0),
      size_(0)
    {}

    partition_data(uint size_x, uint size_y)
    : data_(new T [size_x * size_y], size_x * size_y, buffer_type::take, array_deleter<T>()),
      size_x_(size_x),
      size_y_(size_y),
      size_(size_x * size_y)
    {}

    partition_data(uint size_x, uint size_y, T initial_value)
    : data_(new T [size_x * size_y], size_x * size_y, buffer_type::take, array_deleter<T>()),
      size_x_(size_x),
      size_y_(size_y),
      size_(size_x * size_y)
    {
        for(uint i = 0; i < size_; ++i)
            data_[i] = T(initial_value);
    }

    partition_data(partition_data<T> const& base, partition_type type)
    {
        //return only needed data, depending on who asks for it.
        switch (type)
        {
            case top_left_partition:
            {
                data_ = buffer_type(base.data_.data()+base.size_x()-1, 1, buffer_type::reference);
                size_x_ = 1;
                size_y_ = 1;
                size_ = 1;
                break;
            }

            case top_partition:
            {
                data_ = buffer_type(base.data_.data(), base.size_x(), buffer_type::reference);
                size_x_ = base.size_x();
                size_y_ = 1;
                size_ = base.size_x();
                break;
            }

            case top_right_partition:
            {
                data_ = buffer_type(base.data_.data(), 1, buffer_type::reference);
                size_x_ = 1;
                size_y_ = 1;
                size_ = 1;
                break;
            }

            /*
            * @todo make it so this does not copy
            */
            case left_partition:
            {
                data_ = buffer_type(new T [base.size_y()], base.size_y(), buffer_type::take, array_deleter<T>());
                for(int i = 0; i < base.size_y(); ++i) {
                    data_[i] = base.get_cell(base.size_x()-1,i);
                }

                size_x_ = 1;
                size_y_ = base.size_y();
                size_ = base.size_y();
                break;
            }

            /*
            * @todo make it so this does not copy
            */
            case right_partition:
            {
                data_ = buffer_type(new T [base.size_y()], base.size_y(), buffer_type::take, array_deleter<T>());
                for(int i = 0; i < base.size_y(); ++i) {
                    data_[i] = base.get_cell(0,i);
                }

                size_x_ = 1;
                size_y_ = base.size_y();
                size_ = base.size_y();
                break;
            }

            case bottom_left_partition:
            {
                data_ = buffer_type(base.data_.data()+base.size()-1, 1, buffer_type::reference);
                size_x_ = 1;
                size_y_ = 1;
                size_ = 1;
                break;
            }

            case bottom_partition:
            {
                data_ = buffer_type(base.data_.data()+base.size()-base.size_x(), base.size_x(), buffer_type::reference);
                size_x_ = base.size_x();
                size_y_ = 1;
                size_ = base.size_x();
                break;
            }

            case bottom_right_partition:
            {
                data_ = buffer_type(base.data_.data()+base.size()-base.size_x(), 1, buffer_type::reference);
                size_x_ = 1;
                size_y_ = 1;
                size_ = 1;
                break;
            }

            default:
                HPX_ASSERT(false);
                break;
        }
    }

    uint size_x() const { return size_x_;}
    uint size_y() const { return size_y_;}
    uint size() const { return size_;}

    T get_cell(uint idx, uint idy) const { return data_[index(idx, idy)];}
    T& get_cell_ref(uint idx, uint idy) { return data_[index(idx, idy)];}

    T operator[](uint idx) const { return data_[idx];}
    T& operator[](uint idx) { return data_[idx];}

    friend std::ostream& operator<<(std::ostream& os, partition_data const& data)
    {
        os << "[[";
        for(uint j = 0; j < data.size_y(); ++j)
        {
            for(uint i = 0; i < data.size_x(); ++i)
            {
                os << data.get_cell(i, j).p;
                if(i != data.size_x()-1)
                    os << ", ";
            }
            os << "]";
            if(j < data.size_y()-1)
                os << ",[";
        }
        os << "]";
        return os;
    }

private:
    // Serialization support: even if all of the code below runs on one
    // locality only, we need to provide an (empty) implementation for the
    // serialization as all arguments passed to actions have to support this.
    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & data_ & size_x_ & size_y_ & size_;
    }

    //for accessing an element conveniently
    uint index(uint idx, uint idy) const
    {
        uint id = idy*size_x_ + idx;
        HPX_ASSERT(id >= 0 && id < size_);
        return id;
    }

    buffer_type data_;
    uint size_x_;
    uint size_y_;
    uint size_;
};

}//namespace grid
#endif
