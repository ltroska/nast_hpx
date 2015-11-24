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

    cell(RealType a) { p = a;}

    RealType p;
    RealType u;
    RealType v;
    RealType F;
    RealType G;
};

}

#endif
