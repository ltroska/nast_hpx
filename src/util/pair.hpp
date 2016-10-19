#ifndef NAST_HPX_UTIL_PAIR_HPP_
#define NAST_HPX_UTIL_PAIR_HPP_

namespace nast_hpx {

template<typename T = Real>
struct pair
{
    pair() : x(0), y(0) {}
    pair(T value) : x(value), y(value) {}
    pair(T x_value, T y_value) : x(x_value), y(y_value) {}

    T x, y;

    friend std::ostream& operator<<(std::ostream& os, pair const& p)
    {
        os << "{" << p.x << "," << p.y << "}";

        return os;
    }

    template<typename Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & x & y;
    }

};

typedef pair<std::size_t> index;

}

#endif // PAIR_HPP

