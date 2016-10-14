#ifndef NAST_HPX_GRID_BOUNDARY_DATA_HPP_
#define NAST_HPX_GRID_BOUNDARY_DATA_HPP_

#include "util/defines.hpp"
#include "util/triple.hpp"

#include <iostream>

namespace nast_hpx { namespace grid {


struct data
{

  triple<double> left, right, bottom, top, front, back;
  std::size_t left_type, right_type, bottom_type, top_type, front_type, back_type;

};


/// This class represents one set of boundary conditions.
struct boundary_condition
{
    boundary_condition()
    : left(0), right(0), bottom(0), top(0), back(0), front(0),
      left_type(noslip), right_type(noslip), bottom_type(noslip), top_type(noslip),
      back_type(noslip), front_type(noslip)
    {}

    triple<double> left, right, bottom, top, back, front;
    std::size_t left_type, right_type, bottom_type, top_type, back_type, front_type;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & left & right & bottom & top & back & front
         & left_type & right_type & bottom_type & top_type & back_type & front_type;
    }

    friend std::ostream& operator<<(std::ostream& os, boundary_condition const& data)
    {
        os << "{data: "
                       << "left:" << data.left
                << "," << "right:" << data.right
                << "," << "bottom:" << data.bottom
                << "," << "top:" << data.top
                << "," << "back:" << data.back
                << "," << "front:" << data.front
                << " | types: "
                       << "left:" << data.left_type
                << "," << "right:" << data.right_type
                << "," << "bottom:" << data.bottom_type
                << "," << "top:" << data.top_type
                << "," << "back:" << data.back_type
                << "," << "front:" << data.front_type
                << "}";
        return os;
    }
};

}
}

#endif
