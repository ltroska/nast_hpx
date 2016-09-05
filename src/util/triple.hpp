#ifndef NAST_HPX_UTIL_TRIPLE_HPP
#define NAST_HPX_UTIL_TRIPLE_HPP

namespace nast_hpx {

template<typename T>
struct triple {

    triple() : x(), y(), z() {}

    triple(T value) : x(value), y(value), z(value) {}

    triple(T x_value, T y_value, T z_value) : x(x_value), y(y_value), z(z_value) {}

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
