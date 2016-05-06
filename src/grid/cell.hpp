#ifndef GRID_CELL_HPP
#define GRID_CELL_HPP

#include "util/typedefs.hpp"

namespace grid {

struct scalar_cell
{
    scalar_cell() : value(0) {}

    scalar_cell(RealType v) : value(v) {}

    RealType value;

    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & value;
    }

    friend std::ostream& operator<<(std::ostream& os, scalar_cell const& cell)
    {
        os << "{" << cell  << "}";
        return os;
    }
};

struct vector_cell
{
    vector_cell() : first(0), second(0) {}

    vector_cell(RealType value) : first(value), second(value) {}

    vector_cell(RealType value1, RealType value2)
        : first(value1), second(value2) {}

    RealType first, second;

    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & first & second;
    }

    friend std::ostream& operator<<(std::ostream& os, vector_cell const& cell)
    {
        os << "{" << cell.first << "|" << cell.second << "}";
        return os;
    }
};

}

#endif
