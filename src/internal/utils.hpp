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

inline std::size_t locidx(std::size_t i, std::size_t np, std::size_t nl)
{
    return i / (np/nl);
}


#endif
