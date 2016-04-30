#ifndef GRID_PARTITION_DATA_HPP
#define GRID_PARTITION_DATA_HPP

#include <hpx/runtime/serialization/serialize.hpp>

#include "cell.hpp"
#include "direction.hpp"


namespace grid {

//partition_data
template<typename T = RealType>
struct partition_data
{
private:
    typedef hpx::serialization::serialize_buffer<T> buffer_type;

    struct array_deleter {
        void operator()(T const* p)
        {
            delete [] p;
        }
    };

    struct hold_reference
    {
        hold_reference(buffer_type const& data)
          : data_(data)
        {}

        void operator()(T*) {}     // no deletion necessary

        buffer_type data_;
    };

public:

    partition_data()
    : size_x_(0),
      size_y_(0),
      size_(0)
    {}

    partition_data(uint size_x, uint size_y)
    : data_(new T[size_x*size_y], size_x * size_y, buffer_type::take,
        array_deleter()),
      size_x_(size_x),
      size_y_(size_y),
      size_(size_x * size_y)
    {}

    partition_data(uint size_x, uint size_y, RealType initial_value)
    : data_(new T[size_x*size_y], size_x * size_y, buffer_type::take,
        array_deleter()),
      size_x_(size_x),
      size_y_(size_y),
      size_(size_x * size_y)
    {
        for(uint i = 0; i < size_; ++i)
            data_[i] = T(initial_value);

    }

    partition_data(partition_data const& base, direction type)
    {
        if (base.size() == 0)
        {
            size_x_ = 0;
            size_y_ = 0;
            size_ = 0;
        }
        else
        //return only needed data, depending on who asks for it.
        switch (type)
        {
            case TOP_LEFT:
            {
                data_ =
                    buffer_type(base.data_.data() + base.size_x() - 1, 1,
                                    buffer_type::reference,
                                    hold_reference(base.data_));
                size_x_ = 1;
                size_y_ = 1;
                size_ = 1;
                break;
            }

            case TOP:
            {
                data_ =
                    buffer_type(base.data_.data(), base.size_x(),
                                    buffer_type::reference,
                                    hold_reference(base.data_));
                size_x_ = base.size_x();
                size_y_ = 1;
                size_ = base.size_x();
                break;
            }

            case TOP_RIGHT:
            {
                data_ =
                    buffer_type(base.data_.data(), 1,
                                    buffer_type::reference,
                                    hold_reference(base.data_));
                size_x_ = 1;
                size_y_ = 1;
                size_ = 1;
                break;
            }

            // TODO: make this not copy
            case LEFT:
            {
                data_ =
                    buffer_type(new T [base.size_y()], base.size_y(),
                                    buffer_type::take, array_deleter());
                
                for(uint i = 0; i < base.size_y(); ++i) {
                    data_[i] = base.get_cell(base.size_x() - 1,i);
                }

                size_x_ = 1;
                size_y_ = base.size_y();
                size_ = base.size_y();
                break;
            }

            case RIGHT:
            {
                data_ =
                    buffer_type(new T [base.size_y()], base.size_y(),
                                    buffer_type::take, array_deleter());
                
                for(uint i = 0; i < base.size_y(); ++i) {
                    data_[i] = base.get_cell(0, i);
                }

                size_x_ = 1;
                size_y_ = base.size_y();
                size_ = base.size_y();
                break;
            }

            case BOTTOM_LEFT:
            {
                data_ =
                    buffer_type(base.data_.data() + base.size() - 1, 1,
                                    buffer_type::reference,
                                    hold_reference(base.data_));
                size_x_ = 1;
                size_y_ = 1;
                size_ = 1;
                break;
            }

            case BOTTOM:
            {
                data_ =
                    buffer_type(base.data_.data() + base.size() - base.size_x(),
                                    base.size_x(), buffer_type::reference,
                                    hold_reference(base.data_));
                size_x_ = base.size_x();
                size_y_ = 1;
                size_ = base.size_x();
                break;
            }

            case BOTTOM_RIGHT:
            {
                data_ =
                    buffer_type(base.data_.data() + base.size() - base.size_x(),
                                    1, buffer_type::reference,
                                    hold_reference(base.data_));
                size_x_ = 1;
                size_y_ = 1;
                size_ = 1;
                break;
            }

            case CENTER:
            default:
                data_ =
                    buffer_type(base.data_.data(), base.size(),
                    buffer_type::reference, hold_reference(base.data_));
                size_x_ = base.size_x();
                size_y_ = base.size_y();
                size_ = base.size();
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

    T* begin() { return data_.begin(); }
    T* end() { return data_.end(); }

private:
    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & data_ & size_x_ & size_y_ & size_;
    }

    friend std::ostream& operator<<(std::ostream& os,
        partition_data<T> const& data)
    {
        for (uint j = data.size_y() - 1; j < data.size_y(); j--)
        {
            for (uint i = 0; i < data.size_x(); i++)
                os << data.get_cell(i, j) << " ";

            os << "\n";
        }

        return os;
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
