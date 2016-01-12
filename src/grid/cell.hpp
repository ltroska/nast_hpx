#ifndef GRID_CELL_HPP
#define GRID_CELL_HPP

namespace grid {

struct cell
{
    cell() : p(0), u(0), v(0), f(0), g(0), rhs(0) {}

    cell(RealType value) : p(value), u(value), v(value), f(value), g(value), rhs(value) {}

    RealType p, u, v, f, g, rhs;

    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & p & u & v & f & g & rhs;
    }

    friend std::ostream& operator<<(std::ostream& os, cell const& cell)
    {
        os << "{" << cell.p << "}";
        return os;
    }
};

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
    vector_cell() : u(0), v(0) {}

    vector_cell(RealType value1, RealType value2) : u(value1), v(value2) {}

    RealType u, v;

    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & u & v;
    }
};
}

#endif
