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
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t y, std::size_t z,
                            std::size_t cells_y_per_block, std::size_t cells_z_per_block)
        {
            buffer = BufferType(new Real[cells_y_per_block * cells_z_per_block],
                                    cells_y_per_block * cells_z_per_block, BufferType::take, util::array_deleter<Real>());
            typename BufferType::value_type* src = buffer.data();

            std::size_t start_k = (z == 0 ? 1 : 1 + z * cells_z_per_block);
            std::size_t end_k = start_k + cells_z_per_block;
            std::size_t start_j = (y == 0 ? 1 : 1 + y * cells_y_per_block);
            std::size_t end_j = start_j + cells_y_per_block;

            for (std::size_t k = start_k; k < end_k ; ++k)
                for (std::size_t j = start_j; j < end_j ; ++j)
                {
                    *src = p(1, j, k);
                    ++src;
                }
        }
    };

    template <>
    struct pack_buffer<RIGHT>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t y, std::size_t z,
                            std::size_t cells_y_per_block, std::size_t cells_z_per_block)
        {
            buffer = BufferType(new Real[cells_y_per_block * cells_z_per_block],
                                    cells_y_per_block * cells_z_per_block, BufferType::take, util::array_deleter<Real>());

            typename BufferType::value_type* src = buffer.data();

            std::size_t start_k = (z == 0 ? 1 : 1 + z * cells_z_per_block);
            std::size_t end_k = start_k + cells_z_per_block;
            std::size_t start_j = (y == 0 ? 1 : 1 + y * cells_y_per_block);
            std::size_t end_j = start_j + cells_y_per_block;

            for (std::size_t k = start_k; k < end_k ; ++k)
                for (std::size_t j = start_j; j < end_j ; ++j)
                {
                    *src = p(p.size_x_ - 2, j, k);
                    ++src;
                }
        }
    };

    template <>
    struct pack_buffer<BOTTOM>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t x, std::size_t y,
                            std::size_t cells_x_per_block, std::size_t cells_y_per_block)
        {
            buffer = BufferType(new Real[cells_x_per_block * cells_y_per_block],
                                    cells_x_per_block * cells_y_per_block, BufferType::take, util::array_deleter<Real>());

            typename BufferType::value_type* src = buffer.data();

            std::size_t start_i = (x == 0 ? 1 : 1 + x * cells_x_per_block);
            std::size_t end_i = start_i + cells_x_per_block;
            std::size_t start_j = (y == 0 ? 1 : 1 + y * cells_y_per_block);
            std::size_t end_j = start_j + cells_y_per_block;

            for (std::size_t i = start_i; i < end_i ; ++i)
                for (std::size_t j = start_j; j < end_j ; ++j)
                {
                    *src = p(i, j, 1);
                    ++src;
                }
        }
    };

    template <>
    struct pack_buffer<TOP>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t x, std::size_t y,
                            std::size_t cells_x_per_block, std::size_t cells_y_per_block)
        {
            buffer = BufferType(new Real[cells_x_per_block * cells_y_per_block],
                                    cells_x_per_block * cells_y_per_block, BufferType::take, util::array_deleter<Real>());

            typename BufferType::value_type* src = buffer.data();

            std::size_t start_i = (x == 0 ? 1 : 1 + x * cells_x_per_block);
            std::size_t end_i = start_i + cells_x_per_block;
            std::size_t start_j = (y == 0 ? 1 : 1 + y * cells_y_per_block);
            std::size_t end_j = start_j + cells_y_per_block;

            for (std::size_t i = start_i; i < end_i ; ++i)
                for (std::size_t j = start_j; j < end_j ; ++j)
                {
                    *src = p(i, j, p.size_z_ - 2);
                    ++src;
                }
        }
    };

    template <>
    struct pack_buffer<FRONT>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t x, std::size_t z,
                            std::size_t cells_x_per_block, std::size_t cells_z_per_block)
        {
            buffer = BufferType(new Real[cells_x_per_block * cells_z_per_block],
                                    cells_x_per_block * cells_z_per_block, BufferType::take, util::array_deleter<Real>());

            typename BufferType::value_type* src = buffer.data();

            std::size_t start_i = (x == 0 ? 1 : 1 + x * cells_x_per_block);
            std::size_t end_i = start_i + cells_x_per_block;
            std::size_t start_k = (z == 0 ? 1 : 1 + z * cells_z_per_block);
            std::size_t end_k = start_k + cells_z_per_block;

            for (std::size_t i = start_i; i < end_i; ++i)
                for (std::size_t k = start_k; k < end_k; ++k)
                {
                    *src = p(i, 1, k);
                    ++src;
                }
        }
    };

    template <>
    struct pack_buffer<BACK>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t x, std::size_t z,
                            std::size_t cells_x_per_block, std::size_t cells_z_per_block)
        {
            buffer = BufferType(new Real[cells_x_per_block * cells_z_per_block],
                                    cells_x_per_block * cells_z_per_block, BufferType::take, util::array_deleter<Real>());

            typename BufferType::value_type* src = buffer.data();

            std::size_t start_i = (x == 0 ? 1 : 1 + x * cells_x_per_block);
            std::size_t end_i = start_i + cells_x_per_block;
            std::size_t start_k = (z == 0 ? 1 : 1 + z * cells_z_per_block);
            std::size_t end_k = start_k + cells_z_per_block;

            for (std::size_t i = start_i; i < end_i; ++i)
                for (std::size_t k = start_k; k < end_k; ++k)
                {
                    *src = p(i, p.size_y_ - 2, k);
                    ++src;
                }
        }
    };

    template <>
    struct pack_buffer<BACK_LEFT>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t z, std::size_t bogus,
                            std::size_t cells_z_per_block, std::size_t bogus2)
        {
            buffer = BufferType(new Real[cells_z_per_block],
                                    cells_z_per_block, BufferType::take, util::array_deleter<Real>());

            typename BufferType::value_type* src = buffer.data();

            std::size_t start_k = (z == 0 ? 1 : 1 + z * cells_z_per_block);
            std::size_t end_k = start_k + cells_z_per_block;

            for (std::size_t k = start_k; k < end_k; ++k)
            {
                *src = p(1, p.size_y_ - 2, k);
                ++src;
            }
        }
    };

    template <>
    struct pack_buffer<FRONT_RIGHT>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t z, std::size_t bogus,
                            std::size_t cells_z_per_block, std::size_t bogus2)
        {
            buffer = BufferType(new Real[cells_z_per_block],
                                    cells_z_per_block, BufferType::take, util::array_deleter<Real>());

            typename BufferType::value_type* src = buffer.data();

            std::size_t start_k = (z == 0 ? 1 : 1 + z * cells_z_per_block);
            std::size_t end_k = start_k + cells_z_per_block;

            for (std::size_t k = start_k; k < end_k; ++k)
            {
                *src = p(p.size_x_ - 2, 1, k);
                ++src;
            }
        }
    };

    template <>
    struct pack_buffer<BOTTOM_RIGHT>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t y, std::size_t bogus,
                            std::size_t cells_y_per_block, std::size_t bogus2)
        {
            buffer = BufferType(new Real[cells_y_per_block],
                                    cells_y_per_block, BufferType::take, util::array_deleter<Real>());

            typename BufferType::value_type* src = buffer.data();

            std::size_t start_j = (y == 0 ? 1 : 1 + y * cells_y_per_block);
            std::size_t end_j = start_j + cells_y_per_block;

            for (std::size_t j = start_j; j < end_j; ++j)
            {
                *src = p(p.size_x_ - 2, j, 1);
                ++src;
            }
        }
    };

    template <>
    struct pack_buffer<TOP_LEFT>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t y, std::size_t bogus,
                            std::size_t cells_y_per_block, std::size_t bogus2)
        {
            buffer = BufferType(new Real[cells_y_per_block],
                                    cells_y_per_block, BufferType::take, util::array_deleter<Real>());

            typename BufferType::value_type* src = buffer.data();

            std::size_t start_j = (y == 0 ? 1 : 1 + y * cells_y_per_block);
            std::size_t end_j = start_j + cells_y_per_block;

            for (std::size_t j = start_j; j < end_j; ++j)
            {
                *src = p(1, j, p.size_z_ - 2);
                ++src;
            }
        }
    };

    template <>
    struct pack_buffer<BACK_BOTTOM>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t x, std::size_t bogus,
                            std::size_t cells_x_per_block, std::size_t bogus2)
        {
            buffer = BufferType(new Real[cells_x_per_block],
                                    cells_x_per_block, BufferType::take, util::array_deleter<Real>());

            typename BufferType::value_type* src = buffer.data();

            std::size_t start_i = (x == 0 ? 1 : 1 + x * cells_x_per_block);
            std::size_t end_i = start_i + cells_x_per_block;

            for (std::size_t i = start_i; i < end_i; ++i)
            {
                *src = p(i, p.size_y_ - 2, 1);
                ++src;
            }
        }
    };

    template <>
    struct pack_buffer<FRONT_TOP>
    {
        template <typename BufferType>
        static void call(partition_data<Real> const& p, BufferType& buffer, std::size_t x, std::size_t bogus,
                            std::size_t cells_x_per_block, std::size_t bogus2)
        {
            buffer = BufferType(new Real[cells_x_per_block],
                                    cells_x_per_block, BufferType::take, util::array_deleter<Real>());

            typename BufferType::value_type* src = buffer.data();

            std::size_t start_i = (x == 0 ? 1 : 1 + x * cells_x_per_block);
            std::size_t end_i = start_i + cells_x_per_block;

            for (std::size_t i = start_i; i < end_i; ++i)
            {
                *src = p(i, 1, p.size_z_ - 2);
                ++src;
            }
        }
    };


}
}

#endif
