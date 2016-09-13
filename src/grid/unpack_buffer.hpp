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
        static void call(partition_data<Real>& p, BufferType buffer, std::size_t y, std::size_t z,
                            std::size_t cells_y_per_block, std::size_t cells_z_per_block)
        {
            typename BufferType::value_type* src = buffer.data();

            //HPX_ASSERT(buffer.size() == p.size_y_ - 2);

            std::size_t start_k = (z == 0 ? 1 : 1 + z * cells_z_per_block);
            std::size_t end_k = start_k + cells_z_per_block;
            std::size_t start_j = (y == 0 ? 1 : 1 + y * cells_y_per_block);
            std::size_t end_j = start_j + cells_y_per_block;

            for (std::size_t k = start_k; k < end_k ; ++k)
                for (std::size_t j = start_j; j < end_j ; ++j)
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
        static void call(partition_data<Real>& p, BufferType buffer, std::size_t y, std::size_t z,
                            std::size_t cells_y_per_block, std::size_t cells_z_per_block)
        {
            typename BufferType::value_type* src = buffer.data();

            //HPX_ASSERT(buffer.size() == p.size_y_ - 2);

            std::size_t start_k = (z == 0 ? 1 : 1 + z * cells_z_per_block);
            std::size_t end_k = start_k + cells_z_per_block;
            std::size_t start_j = (y == 0 ? 1 : 1 + y * cells_y_per_block);
            std::size_t end_j = start_j + cells_y_per_block;

            for (std::size_t k = start_k; k < end_k ; ++k)
                for (std::size_t j = start_j; j < end_j ; ++j)
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
        static void call(partition_data<Real>& p, BufferType& buffer, std::size_t x, std::size_t y,
                            std::size_t cells_x_per_block, std::size_t cells_y_per_block)
        {
            typename BufferType::value_type* src = buffer.data();

            std::size_t start_i = (x == 0 ? 1 : 1 + x * cells_x_per_block);
            std::size_t end_i = start_i + cells_x_per_block;
            std::size_t start_j = (y == 0 ? 1 : 1 + y * cells_y_per_block);
            std::size_t end_j = start_j + cells_y_per_block;

            for (std::size_t i = start_i; i < end_i ; ++i)
                for (std::size_t j = start_j; j < end_j ; ++j)
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
        static void call(partition_data<Real>& p, BufferType& buffer, std::size_t x, std::size_t y,
                            std::size_t cells_x_per_block, std::size_t cells_y_per_block)
        {
            typename BufferType::value_type* src = buffer.data();

            std::size_t start_i = (x == 0 ? 1 : 1 + x * cells_x_per_block);
            std::size_t end_i = start_i + cells_x_per_block;
            std::size_t start_j = (y == 0 ? 1 : 1 + y * cells_y_per_block);
            std::size_t end_j = start_j + cells_y_per_block;

            for (std::size_t i = start_i; i < end_i; ++i)
                for (std::size_t j = start_j; j < end_j; ++j)
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
        static void call(partition_data<Real>& p, BufferType& buffer, std::size_t x, std::size_t z,
                            std::size_t cells_x_per_block, std::size_t cells_z_per_block)
        {
            typename BufferType::value_type* src = buffer.data();

            std::size_t start_i = (x == 0 ? 1 : 1 + x * cells_x_per_block);
            std::size_t end_i = start_i + cells_x_per_block;
            std::size_t start_k = (z == 0 ? 1 : 1 + z * cells_z_per_block);
            std::size_t end_k = start_k + cells_z_per_block;

            for (std::size_t i = start_i; i < end_i; ++i)
                for (std::size_t k = start_k; k < end_k; ++k)
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
        static void call(partition_data<Real>& p, BufferType& buffer, std::size_t x, std::size_t z,
                            std::size_t cells_x_per_block, std::size_t cells_z_per_block)
        {
            typename BufferType::value_type* src = buffer.data();

            std::size_t start_i = (x == 0 ? 1 : 1 + x * cells_x_per_block);
            std::size_t end_i = start_i + cells_x_per_block;
            std::size_t start_k = (z == 0 ? 1 : 1 + z * cells_z_per_block);
            std::size_t end_k = start_k + cells_z_per_block;

            for (std::size_t i = start_i; i < end_i; ++i)
                for (std::size_t k = start_k; k < end_k; ++k)
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
        static void call(partition_data<Real>& p, BufferType& buffer, std::size_t z, std::size_t bogus,
                            std::size_t cells_z_per_block, std::size_t bogus2)
        {
            typename BufferType::value_type* src = buffer.data();

            std::size_t start_k = (z == 0 ? 1 : 1 + z * cells_z_per_block);
            std::size_t end_k = start_k + cells_z_per_block;

            for (std::size_t k = start_k; k < end_k; ++k)
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
        static void call(partition_data<Real>& p, BufferType& buffer, std::size_t z, std::size_t bogus,
                            std::size_t cells_z_per_block, std::size_t bogus2)
        {
            typename BufferType::value_type* src = buffer.data();

            std::size_t start_k = (z == 0 ? 1 : 1 + z * cells_z_per_block);
            std::size_t end_k = start_k + cells_z_per_block;

            for (std::size_t k = start_k; k < end_k; ++k)
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
        static void call(partition_data<Real>& p, BufferType& buffer, std::size_t y, std::size_t bogus,
                            std::size_t cells_y_per_block, std::size_t bogus2)
        {
            typename BufferType::value_type* src = buffer.data();

            std::size_t start_j = (y == 0 ? 1 : 1 + y * cells_y_per_block);
            std::size_t end_j = start_j + cells_y_per_block;

            for (std::size_t j = start_j; j < end_j; ++j)
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
        static void call(partition_data<Real>& p, BufferType& buffer, std::size_t y, std::size_t bogus,
                            std::size_t cells_y_per_block, std::size_t bogus2)
        {
            typename BufferType::value_type* src = buffer.data();

            std::size_t start_j = (y == 0 ? 1 : 1 + y * cells_y_per_block);
            std::size_t end_j = start_j + cells_y_per_block;

            for (std::size_t j = start_j; j < end_j; ++j)
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
        static void call(partition_data<Real>& p, BufferType& buffer, std::size_t x, std::size_t bogus,
                            std::size_t cells_x_per_block, std::size_t bogus2)
        {
            typename BufferType::value_type* src = buffer.data();

            std::size_t start_i = (x == 0 ? 1 : 1 + x * cells_x_per_block);
            std::size_t end_i = start_i + cells_x_per_block;

            for (std::size_t i = start_i; i < end_i; ++i)
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
        static void call(partition_data<Real>& p, BufferType& buffer, std::size_t x, std::size_t bogus,
                            std::size_t cells_x_per_block, std::size_t bogus2)
        {
            typename BufferType::value_type* src = buffer.data();

            std::size_t start_i = (x == 0 ? 1 : 1 + x * cells_x_per_block);
            std::size_t end_i = start_i + cells_x_per_block;

            for (std::size_t i = start_i; i < end_i; ++i)
            {
                p(i, 0, p.size_z_ - 1) = *src;
                ++src;
            }
        }
    };

}
}

#endif
