#ifndef UTIL_TYPES_HPP
#define UTIL_TYPES_HPP

#define RealType double
#define uint std::size_t

//directions for neighbor relations
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
