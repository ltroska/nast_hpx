#ifndef NAST_HPX_GRID_PACK_BUFFER_HPP
#define NAST_HPX_GRID_PACK_BUFFER_HPP

#include <hpx/runtime/serialization/serialize_buffer.hpp>

#include "direction.hpp"
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

        }
    };

    template <>
    struct pack_buffer<RIGHT>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t offset, std::size_t size)
        {

        }
    };

    template <>
    struct pack_buffer<BOTTOM>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t offset, std::size_t size)
        {
        }
    };

    template <>
    struct pack_buffer<TOP>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t offset, std::size_t size)
        {
        }
    };

    template <>
    struct pack_buffer<BOTTOM_RIGHT>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t offset, std::size_t size)
        {
        }
    };

    template <>
    struct pack_buffer<TOP_LEFT>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t offset, std::size_t size)
        {
        }
    };

}
}

#endif
