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
                        
            for(std::size_t y = 1; y != p.size_y_ - 1; ++y)
            {
                *src = p(1, y);
                ++src;
            }
        }
    };    
    
    template <>
    struct pack_buffer<RIGHT>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer)
        {
            Real* data = new Real[p.act_size_y_];
            buffer = BufferType(data, p.act_size_y_, BufferType::take);
            
            typename BufferType::value_type * src = buffer.data();
                        
            for(std::size_t y = 1; y != p.size_y_ - 1; ++y)
            {
                *src = p(p.size_x_ - 2, y);
                ++src;
            }
        }
    }; 
   
    template <>
    struct pack_buffer<BOTTOM>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer)
        {
            buffer = BufferType(p.data_.data() + p.size_x_ + 1, p.act_size_x_, BufferType::reference);
        }
    };    
    
    template <>
    struct pack_buffer<TOP>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer)
        {
            buffer = BufferType(p.data_.data() + p.size_ - 2 * p.size_x_ + 1, p.act_size_x_, BufferType::reference);
        }
    };    
        
    template <>
    struct pack_buffer<BOTTOM_RIGHT>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer)
        {
            buffer = BufferType(p.data_.data() + 2 * p.size_x_ - 2, 1, BufferType::reference);
        }
    };
    
    template <>
    struct pack_buffer<TOP_LEFT>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer)
        {
            buffer = BufferType(p.data_.data() + p.size_ - 2 * p.size_x_ + 1, 1, BufferType::reference);
        }
    };    

}
}

#endif
