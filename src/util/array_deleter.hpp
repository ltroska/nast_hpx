#ifndef NAST_HPX_UTIL_ARRAY_DELETER_HPP_
#define NAST_HPX_UTIL_ARRAY_DELETER_HPP_

namespace nast_hpx { namespace util {

template<typename T>
struct array_deleter {
    void operator()(T const* p)
    {
        delete [] p;
    }
};

}
}

#endif
