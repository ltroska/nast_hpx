#ifndef GRID_CELL_HPP
#define GRID_CELL_HPP

namespace grid {

struct cell {
public:
    cell()
    : p(0),
      u(0),
      v(0),
      F(0),
      G(0)
    {}

    cell(RealType a)
    {
        cell();
        p = a;
    }

    RealType p;
    RealType u;
    RealType v;
    RealType F;
    RealType G;

    friend std::ostream& operator<<(std::ostream& os, cell const& cell)
    {
        os << "{" << cell.p << "|" << cell.u << "|" << cell.v << "}";
        return os;
    }

    friend class hpx::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & p & u & v & F & G;
    }
};

}

#endif
