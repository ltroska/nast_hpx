#ifndef NAST_HPX_UTIL_TRIPLE_HPP_
#define NAST_HPX_UTIL_TRIPLE_HPP_

#include <iostream>

namespace nast_hpx {

template<typename T>
struct triple {

    triple() : x(), y(), z() {}

    triple(T value) : x(value), y(value), z(value) {}

    triple(T x_value, T y_value, T z_value) : x(x_value), y(y_value), z(z_value) {}

    triple(triple<T> const& other) : x(other.x), y(other.y), z(other.z) {}

    triple(triple<T>&& other) : x(std::move(other.x)), y(std::move(other.y)), z(std::move(other.z)) {}

    triple<T>& operator=(triple<T> const& other)
    {
        x = other.x;
        y = other.y;
        z = other.z;

        return *this;
    }

    triple<T>& operator=(triple<T>&& other)
    {
        x = std::move(other.x);
        y = std::move(other.y);
        z = std::move(other.z);

        return *this;
    }

    T x, y, z;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & x & y & z;
    }

    friend std::ostream& operator<<(std::ostream& os, triple<T> const& triple)
    {
        os << "{" << triple.x << "," << triple.y << "," << triple.z << "}";
        return os;
    }

};

typedef triple<std::size_t> index;

}

#endif
