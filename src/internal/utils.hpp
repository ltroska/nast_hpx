#ifndef INTERNAL_UTILS_HPP
#define INTERNAL_UTILS_HPP

//needed for deleting the cell arrays properly
template<typename T>
class array_deleter {
    public:
        void operator()(T const* p)
        {
            delete [] p;
        }
};

#endif
