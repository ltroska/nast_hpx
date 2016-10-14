#ifndef NAST_HPX_GRID_DIRECTION_HPP_
#define NAST_HPX_GRID_DIRECTION_HPP_

namespace nast_hpx { namespace grid {

/// This enum represents all directions needed for 3D neighbor relations
enum direction
{
    LEFT = 0,
    BOTTOM,
    BACK,
    BACK_LEFT,
    BOTTOM_RIGHT,
    BACK_BOTTOM,
    FRONT_TOP,
    TOP_LEFT,
    FRONT_RIGHT,
    FRONT,
    TOP,
    RIGHT,
    NUM_DIRECTIONS
};

}
}
#endif
