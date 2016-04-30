#pragma once
#ifndef BOUNDARY_DATA_HPP_
#define BOUNDARY_DATA_HPP_

#include <iostream>

/// This class represents one set of boundary conditions.
struct boundary_data
{
    boundary_data() : left(0), right(0), bottom(0), top(0) {}

    boundary_data(RealType value)
    : left(value), right(value), bottom(value), top(value) {}

    boundary_data(RealType left_value, RealType right_value,
        RealType bottom_value, RealType top_value)
    : left(left_value), right(right_value), bottom(bottom_value),
        top(top_value) {}

    RealType left, right, bottom, top;

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

#endif
