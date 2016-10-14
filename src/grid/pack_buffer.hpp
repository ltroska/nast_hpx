#ifndef NAST_HPX_GRID_PACK_BUFFER_HPP_
#define NAST_HPX_GRID_PACK_BUFFER_HPP_

#include "direction.hpp"
#include "partition_data.hpp"
#include "util/array_deleter.hpp"

#include <hpx/runtime/serialization/serialize_buffer.hpp>

namespace nast_hpx { namespace grid {
    template <direction dir>
    struct pack_buffer;

    template <>
    struct pack_buffer<LEFT>
    {
        template <typename BufferType>
        static void call(partition_data<double> const& p, BufferType& buffer)
        {
            buffer = BufferType(new double[(p.size_z_ - 2) * (p.size_y_ - 2)],
                                    (p.size_z_ - 2) * (p.size_y_ - 2), BufferType::take, util::array_deleter<double>());
            typename BufferType::value_type* src = buffer.data();


            for (std::size_t k = 1; k < p.size_z_ - 1; ++k)
                for (std::size_t j = 1; j < p.size_y_ - 1; ++j)
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
        static void call(partition_data<double> const& p, BufferType& buffer)
        {
            buffer = BufferType(new double[(p.size_z_ - 2) * (p.size_y_ - 2)],
                                   (p.size_z_ - 2) * (p.size_y_ - 2), BufferType::take, util::array_deleter<double>());

            typename BufferType::value_type* src = buffer.data();

            for (std::size_t k = 1; k < p.size_z_ - 1; ++k)
                for (std::size_t j = 1; j < p.size_y_ - 1; ++j)
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
        static void call(partition_data<double> const& p, BufferType& buffer)
        {
            buffer = BufferType(new double[(p.size_x_ - 2) * (p.size_y_ - 2)],
                                    (p.size_x_ - 2) * (p.size_y_ - 2), BufferType::take, util::array_deleter<double>());

            typename BufferType::value_type* src = buffer.data();

            for (std::size_t i = 1; i < p.size_x_ - 1; ++i)
                for (std::size_t j = 1; j < p.size_y_ - 1; ++j)
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
        static void call(partition_data<double> const& p, BufferType& buffer)
        {
            buffer = BufferType(new double[(p.size_x_ - 2) * (p.size_y_ - 2)],
                                    (p.size_x_ - 2) * (p.size_y_ - 2), BufferType::take, util::array_deleter<double>());

            typename BufferType::value_type* src = buffer.data();

            for (std::size_t i = 1; i < p.size_x_ - 1; ++i)
                for (std::size_t j = 1; j < p.size_y_ - 1; ++j)
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
        static void call(partition_data<double> const& p, BufferType& buffer)
        {
            buffer = BufferType(new double[(p.size_x_ - 2) * (p.size_z_ - 2)],
                                    (p.size_x_ - 2) * (p.size_z_ - 2), BufferType::take, util::array_deleter<double>());

            typename BufferType::value_type* src = buffer.data();

            for (std::size_t i = 1; i < p.size_x_ - 1; ++i)
                for (std::size_t k = 1; k < p.size_z_ - 1; ++k)
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
        static void call(partition_data<double> const& p, BufferType& buffer)
        {
            buffer = BufferType(new double[(p.size_x_ - 2) * (p.size_z_ - 2)],
                                   (p.size_x_ - 2) * (p.size_z_ - 2), BufferType::take, util::array_deleter<double>());

            typename BufferType::value_type* src = buffer.data();

            for (std::size_t i = 1; i < p.size_x_ - 1; ++i)
                for (std::size_t k = 1; k < p.size_z_ - 1; ++k)
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
        static void call(partition_data<double> const& p, BufferType& buffer)
        {
            buffer = BufferType(new double[p.size_z_ - 2],
                                    p.size_z_ - 2, BufferType::take, util::array_deleter<double>());

            typename BufferType::value_type* src = buffer.data();

            for (std::size_t k = 1; k < p.size_z_ - 1; ++k)
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
        static void call(partition_data<double> const& p, BufferType& buffer)
        {
            buffer = BufferType(new double[p.size_z_ - 2],
                                    p.size_z_ - 2, BufferType::take, util::array_deleter<double>());

            typename BufferType::value_type* src = buffer.data();

            for (std::size_t k = 1; k < p.size_z_ - 1; ++k)
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
        static void call(partition_data<double> const& p, BufferType& buffer)
        {
            buffer = BufferType(new double[p.size_y_ - 2],
                                    p.size_y_ - 2, BufferType::take, util::array_deleter<double>());

            typename BufferType::value_type* src = buffer.data();

            for (std::size_t j = 1; j < p.size_y_ - 1; ++j)
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
        static void call(partition_data<double> const& p, BufferType& buffer)
        {
            buffer = BufferType(new double[p.size_y_ - 2],
                                    p.size_y_ - 2, BufferType::take, util::array_deleter<double>());

            typename BufferType::value_type* src = buffer.data();

            for (std::size_t j = 1; j < p.size_y_ - 1; ++j)
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
        static void call(partition_data<double> const& p, BufferType& buffer)
        {
            buffer = BufferType(new double[p.size_x_ - 2],
                                    p.size_x_ - 2, BufferType::take, util::array_deleter<double>());

            typename BufferType::value_type* src = buffer.data();

            for (std::size_t i = 1; i < p.size_x_ - 1; ++i)
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
        static void call(partition_data<double> const& p, BufferType& buffer)
        {
            buffer = BufferType(new double[p.size_x_ - 2],
                                    p.size_x_ - 2, BufferType::take, util::array_deleter<double>());

            typename BufferType::value_type* src = buffer.data();

            for (std::size_t i = 1; i < p.size_z_ - 1; ++i)
            {
                *src = p(i, 1, p.size_z_ - 2);
                ++src;
            }
        }
    };


}
}

#endif
