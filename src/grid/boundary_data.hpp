#ifndef NAST_HPX_GRID_BOUNDARY_DATA_HPP_
#define NAST_HPX_GRID_BOUNDARY_DATA_HPP_

#include <iostream>

#include "util/typedefs.hpp"

namespace nast_hpx { namespace grid {

/// This class represents one set of boundary conditions.
struct boundary_data
{
    boundary_data() : left(0), right(0), bottom(0), top(0) {}

    boundary_data(Real value)
    : left(value), right(value), bottom(value), top(value) {}

    boundary_data(Real left_value, Real right_value,
        Real bottom_value, Real top_value)
    : left(left_value), right(right_value), bottom(bottom_value),
        top(top_value) {}

    Real left, right, bottom, top;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & left & right & bottom & top;
    }

    friend std::ostream& operator<<(std::ostream& os, boundary_data const& data)
    {
        os << "boundary_data={left:" << data.left
                << "|" << "right:" << data.right
                << "|" << "bottom:" << data.bottom
                << "|" << "top:" << data.top << "}";
        return os;
    }
};

}
}

#endif
