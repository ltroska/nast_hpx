#ifndef NAST_HPX_GRID_DIRECTION_HPP
#define NAST_HPX_GRID_DIRECTION_HPP

namespace nast_hpx { namespace grid {

/// This enum represents all directions needed for 3D neighbor relations
enum direction
{
    LEFT = 0,
    BOTTOM,
    BACK,
    BOTTOM_LEFT,
    BOTTOM_RIGHT,
    BACK_BOTTOM,
    BACK_TOP,
    BOTTOM_LEFT_BACK,
    BOTTOM_RIGHT_BACK,
    TOP_LEFT_BACK,
    TOP_RIGHT_BACK,
    BOTTOM_LEFT_FRONT,
    BOTTOM_RIGHT_FRONT,
    TOP_LEFT_FRONT,
    TOP_RIGHT_FRONT,
    FRONT_BOTTOM,
    FRONT_TOP,
    TOP_LEFT,
    TOP_RIGHT,
    FRONT,
    TOP,
    RIGHT,
    NUM_DIRECTIONS
};

}
}
#endif
