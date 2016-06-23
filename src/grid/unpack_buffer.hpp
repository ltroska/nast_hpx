#ifndef NAST_HPX_GRID_UNPACK_BUFFER_HPP
#define NAST_HPX_GRID_UNPACK_BUFFER_HPP

#include <hpx/runtime/serialization/serialize_buffer.hpp>

#include "partition_data.hpp"

namespace nast_hpx { namespace grid {
    template <direction dir>
    struct unpack_buffer;

    template <>
    struct unpack_buffer<LEFT>
    {
        template <typename BufferType>
        static void call(partition_data<Real>& p, BufferType buffer, std::size_t offset)
        {
            typename BufferType::value_type* src = buffer.data();

            HPX_ASSERT(buffer.size() == p.size_y_ - 2);


            for(std::size_t y = 1 + offset; y != 1 + offset + buffer.size(); ++y)
            {
                p(0, y) = *src;
                ++src;
            }
        }
    };

    template <>
    struct unpack_buffer<RIGHT>
    {
        template <typename BufferType>
        static void call(partition_data<Real>& p, BufferType buffer, std::size_t offset)
        {
            typename BufferType::value_type* src = buffer.data();

            HPX_ASSERT(buffer.size() == p.size_y_ - 2);


            for(std::size_t y = 1 + offset; y != 1 + offset + buffer.size(); ++y)
            {
                p(p.size_x_ - 1, y) = *src;
                ++src;
            }
        }
    };

    template <>
    struct unpack_buffer<BOTTOM>
    {
        template <typename BufferType>
        static void call(partition_data<Real>& p, BufferType buffer, std::size_t offset)
        {
            typename BufferType::value_type* src = buffer.data();

            HPX_ASSERT(buffer.size() == p.size_x_ - 2);


            for(std::size_t x = 1 + offset; x != 1 + offset + buffer.size(); ++x)
            {
                p(x, 0) = *src;
                ++src;
            }
        }
    };

    template <>
    struct unpack_buffer<TOP>
    {
        template <typename BufferType>
        static void call(partition_data<Real>& p, BufferType buffer, std::size_t offset)
        {
            typename BufferType::value_type* src = buffer.data();

            HPX_ASSERT(buffer.size() == p.size_x_ - 2);


            for(std::size_t x = 1 + offset; x != 1 + offset + buffer.size(); ++x)
            {
                p(x, p.size_y_ - 1) = *src;
                ++src;
            }
        }
    };

    template <>
    struct unpack_buffer<BOTTOM_RIGHT>
    {
        template <typename BufferType>
        static void call(partition_data<Real>& p, BufferType buffer, std::size_t offset)
        {
            typename BufferType::value_type* src = buffer.data();

            HPX_ASSERT(buffer.size() == 1);

            p(p.size_x_ - 1, 0) = *src;
        }
    };

    template <>
    struct unpack_buffer<TOP_LEFT>
    {
        template <typename BufferType>
        static void call(partition_data<Real>& p, BufferType buffer, std::size_t offset)
        {
            typename BufferType::value_type* src = buffer.data();

            HPX_ASSERT(buffer.size() == 1);

            p(0, p.size_y_ - 1) = *src;
        }
    };

}
}

#endif
