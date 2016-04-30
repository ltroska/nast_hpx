#ifndef GRID_DIRECTION_HPP
#define GRID_DIRECTION_HPP

/// This enum represents all directions needed for 2D neighbor relations
enum direction
{
    LEFT = 0,
    TOP,
    BOTTOM_LEFT,
    BOTTOM_RIGHT,
    CENTER,
    TOP_LEFT,
    TOP_RIGHT,
    BOTTOM,
    RIGHT,
    NUM_DIRECTIONS
};

#endif
