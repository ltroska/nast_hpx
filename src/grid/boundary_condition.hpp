#ifndef NAST_HPX_GRID_BOUNDARY_CONDITION_HPP_
#define NAST_HPX_GRID_BOUNDARY_CONDITION_HPP_

#include "../util/defines.hpp"
#include "../util/pair.hpp"

#include <iostream>

namespace nast_hpx { namespace grid {

/// This class represents one set of boundary conditions.
struct boundary_condition
{
    boundary_condition()
    : left(0), right(0), bottom(0), top(0),
      left_type(noslip), right_type(noslip), bottom_type(noslip), top_type(noslip)
    {}

    pair<Real> left, right, bottom, top;
    std::size_t left_type, right_type, bottom_type, top_type;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & left & right & bottom & top
         & left_type & right_type & bottom_type & top_type;
    }

    friend std::ostream& operator<<(std::ostream& os, boundary_condition const& data)
    {
        os << "{data: "
                       << "left:" << data.left
                << "," << "right:" << data.right
                << "," << "bottom:" << data.bottom
                << "," << "top:" << data.top
                << " | types: "
                       << "left:" << data.left_type
                << "," << "right:" << data.right_type
                << "," << "bottom:" << data.bottom_type
                << "," << "top:" << data.top_type
                << "}";
        return os;
    }
};

}
}

#endif
