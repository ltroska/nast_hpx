#ifndef NAST_HPX_GRID_BOUNDARY_DATA_HPP_
#define NAST_HPX_GRID_BOUNDARY_DATA_HPP_

#include <iostream>

#include "util/typedefs.hpp"

namespace nast_hpx { namespace grid {

/// This class represents one set of boundary conditions.
template<typename T>
struct data
{
    data() : left(0), right(0), bottom(0), top(0) {}

    data(T value)
    : left(value), right(value), bottom(value), top(value) {}

    data(T left_value, T right_value,
        T bottom_value, T top_value)
    : left(left_value), right(right_value), bottom(bottom_value),
        top(top_value) {}

    T left, right, bottom, top;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & left & right & bottom & top;
    }

    friend std::ostream& operator<<(std::ostream& os, data const& data)
    {
        os << "data={left:" << data.left
                << "|" << "right:" << data.right
                << "|" << "bottom:" << data.bottom
                << "|" << "top:" << data.top << "}";
        return os;
    }
};

typedef data<Real> boundary_data;
typedef data<std::size_t> boundary_type;

}
}

#endif
