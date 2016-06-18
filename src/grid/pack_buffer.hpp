#ifndef NAST_HPX_GRID_PACK_BUFFER_HPP
#define NAST_HPX_GRID_PACK_BUFFER_HPP

#include <hpx/runtime/serialization/serialize_buffer.hpp>

#include "partition_data.hpp"
#include "util/array_deleter.hpp"

namespace nast_hpx { namespace grid {
    template <direction dir>
    struct pack_buffer;

    template <>
    struct pack_buffer<LEFT>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t offset, std::size_t size)
        {
            buffer = BufferType(new Real[size], size, BufferType::take, util::array_deleter<Real>());

            typename BufferType::value_type * src = buffer.data();

            for(std::size_t y = 1 + offset; y != 1 + offset + size; ++y)
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
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t offset, std::size_t size)
        {
            buffer = BufferType(new Real[size], size, BufferType::take, util::array_deleter<Real>());

            typename BufferType::value_type * src = buffer.data();

            for(std::size_t y = 1 + offset; y != 1 + offset + size; ++y)
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
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t offset, std::size_t size)
        {
            buffer = BufferType(p.data_.data() + p.size_x_ + 1 + offset, size, BufferType::reference);
        }
    };

    template <>
    struct pack_buffer<TOP>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t offset, std::size_t size)
        {
            buffer = BufferType(p.data_.data() + p.size_ - 2 * p.size_x_ + 1 + offset, size, BufferType::reference);
        }
    };

    template <>
    struct pack_buffer<BOTTOM_RIGHT>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t offset, std::size_t size)
        {
            buffer = BufferType(p.data_.data() + 2 * p.size_x_ - 2, 1, BufferType::reference);
        }
    };

    template <>
    struct pack_buffer<TOP_LEFT>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t offset, std::size_t size)
        {
            buffer = BufferType(p.data_.data() + p.size_ - 2 * p.size_x_ + 1, 1, BufferType::reference);
        }
    };

}
}

#endif
