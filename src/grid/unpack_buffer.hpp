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
        static void call(partition_data<Real>& p, BufferType buffer)
        {
            typename BufferType::value_type* src = buffer.data();

            //HPX_ASSERT(buffer.size() == p.size_y_ - 2);

            for (std::size_t k = 1; k < p.size_z_ - 1 ; ++k)
                for (std::size_t j = 1; j < p.size_y_ - 1 ; ++j)
                {
                    p(0, j, k) = *src;
                    ++src;
                }
        }
    };

    template <>
    struct unpack_buffer<RIGHT>
    {
        template <typename BufferType>
        static void call(partition_data<Real>& p, BufferType buffer)
        {
            typename BufferType::value_type* src = buffer.data();

            //HPX_ASSERT(buffer.size() == p.size_y_ - 2);

            for (std::size_t k = 1; k < p.size_z_ - 1 ; ++k)
                for (std::size_t j = 1; j < p.size_y_ - 1 ; ++j)
                {
                    p(p.size_x_ - 1, j, k) = *src;
                    ++src;
                }
        }
    };

    template <>
    struct unpack_buffer<BOTTOM>
    {
        template <typename BufferType>
        static void call(partition_data<Real>& p, BufferType& buffer)
        {
            typename BufferType::value_type* src = buffer.data();

            for (std::size_t i = 1; i < p.size_x_ - 1 ; ++i)
                for (std::size_t j = 1; j < p.size_y_ - 1 ; ++j)
                {
                    p(i, j, 0) = *src;
                    ++src;
                }
        }
    };

    template <>
    struct unpack_buffer<TOP>
    {
        template <typename BufferType>
        static void call(partition_data<Real>& p, BufferType& buffer)
        {
            typename BufferType::value_type* src = buffer.data();

            for (std::size_t i = 1; i < p.size_x_ - 1; ++i)
                for (std::size_t j = 1; j < p.size_y_ - 1; ++j)
                {
                    p(i, j, p.size_z_ - 1) = *src;
                    ++src;
                }
        }
    };

    template <>
    struct unpack_buffer<FRONT>
    {
        template <typename BufferType>
        static void call(partition_data<Real>& p, BufferType& buffer)
        {
            typename BufferType::value_type* src = buffer.data();

            for (std::size_t i = 1; i < p.size_x_ - 1; ++i)
                for (std::size_t k = 1; k < p.size_z_ - 1; ++k)
                {
                    p(i, 0, k) = *src;
                    ++src;
                }
        }
    };

    template <>
    struct unpack_buffer<BACK>
    {
        template <typename BufferType>
        static void call(partition_data<Real>& p, BufferType& buffer)
        {
            typename BufferType::value_type* src = buffer.data();

            for (std::size_t i = 1; i < p.size_x_ - 1; ++i)
                for (std::size_t k = 1; k < p.size_z_ - 1; ++k)
                {
                    p(i, p.size_y_ - 1, k) = *src;
                    ++src;
                }
        }
    };

    template <>
    struct unpack_buffer<BACK_LEFT>
    {
        template <typename BufferType>
        static void call(partition_data<Real>& p, BufferType& buffer)
        {
            typename BufferType::value_type* src = buffer.data();

            for (std::size_t k = 1; k < p.size_z_ - 1; ++k)
            {
                p(0, p.size_y_ - 1, k) = *src;
                ++src;
            }
        }
    };

    template <>
    struct unpack_buffer<FRONT_RIGHT>
    {
        template <typename BufferType>
        static void call(partition_data<Real>& p, BufferType& buffer)
        {
            typename BufferType::value_type* src = buffer.data();

            for (std::size_t k = 1; k < p.size_z_ - 1; ++k)
            {
                p(p.size_x_ - 1, 0, k) = *src;
                ++src;
            }
        }
    };

    template <>
    struct unpack_buffer<BOTTOM_RIGHT>
    {
        template <typename BufferType>
        static void call(partition_data<Real>& p, BufferType& buffer)
        {
            typename BufferType::value_type* src = buffer.data();

            for (std::size_t j = 1; j < p.size_y_ - 1; ++j)
            {
                p(p.size_x_ - 1, j, 0) = *src;
                ++src;
            }
        }
    };
    template <>
    struct unpack_buffer<TOP_LEFT>
    {
        template <typename BufferType>
        static void call(partition_data<Real>& p, BufferType& buffer)
        {
            typename BufferType::value_type* src = buffer.data();

            for (std::size_t j = 1; j < p.size_y_ - 1; ++j)
            {
                p(0, j, p.size_z_ - 1) = *src;
                ++src;
            }
        }
    };

    template <>
    struct unpack_buffer<BACK_BOTTOM>
    {
        template <typename BufferType>
        static void call(partition_data<Real>& p, BufferType& buffer)
        {
            typename BufferType::value_type* src = buffer.data();

            for (std::size_t i = 1; i < p.size_x_ - 1; ++i)
            {
                p(i, p.size_y_ - 1, 0) = *src;
                ++src;
            }
        }
    };

    template <>
    struct unpack_buffer<FRONT_TOP>
    {
        template <typename BufferType>
        static void call(partition_data<Real>& p, BufferType& buffer)
        {
            typename BufferType::value_type* src = buffer.data();

            for (std::size_t i = 1; i < p.size_x_ - 1; ++i)
            {
                p(i, 0, p.size_z_ - 1) = *src;
                ++src;
            }
        }
    };

}
}

#endif
