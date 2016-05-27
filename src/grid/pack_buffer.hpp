#ifndef NAST_HPX_GRID_PACK_BUFFER_HPP
#define NAST_HPX_GRID_PACK_BUFFER_HPP

#include <hpx/runtime/serialization/serialize_buffer.hpp>

#include "partition_data.hpp"

namespace nast_hpx { namespace grid {
    template <direction dir>
    struct pack_buffer;

    template <>
    struct pack_buffer<LEFT>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer)
        {
            Real* data = new Real[p.act_size_y_];
            buffer = BufferType(data, p.act_size_y_, BufferType::take);
            
            typename BufferType::value_type * src = buffer.data();
            
            std::cout << "packing";
            
            for(std::size_t y = 1; y != p.size_y_ - 1; ++y)
            {
                std::cout << p(1, y) << "\n";
                *src = p(1, y);
                ++src;
            }
            
            std::cout << std::endl;
        }
    };    

}
}

#endif
