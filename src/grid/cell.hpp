#ifndef GRID_CELL_HPP
#define GRID_CELL_HPP

namespace grid {

struct scalar_cell
{
    scalar_cell() : c(0) {}

    scalar_cell(RealType value) : c(value) {}

    RealType c;

    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & c;
    }

    friend std::ostream& operator<<(std::ostream& os, scalar_cell const& cell)
    {
        os << "{" << cell.c << "}";
        return os;
    }
};

struct vector_cell
{
    vector_cell() : c1(0), c2(0) {}

    vector_cell(RealType value) : c1(value), c2(value) {}

    vector_cell(RealType value1, RealType value2) : c1(value1), c2(value2) {}

    RealType c1, c2;

    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & c1 & c2;
    }

    friend std::ostream& operator<<(std::ostream& os, vector_cell const& cell)
    {
        os << "{" << cell.c1 << "|" << cell.c2 << "}";
        return os;
    }
};
}

#endif
