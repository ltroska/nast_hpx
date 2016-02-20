#ifndef GRID_SERVER_PARTITION_SERVER_HPP
#define GRID_SERVER_PARTITION_SERVER_HPP

#include <hpx/include/components.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/util.hpp>

#include <boost/preprocessor/cat.hpp>

#include "grid/partition_data.hpp"
#include "util/types.hpp"

namespace grid { namespace server {

template<typename T>
struct HPX_COMPONENT_EXPORT partition_server
    : public hpx::components::locking_hook<
            hpx::components::simple_component_base<partition_server<T> > >
{
    partition_server() {}

    partition_server(partition_data<T> const& data)
    : data_(data)
    {}

    partition_server(uint size_x, uint size_y, RealType initial_value)
    : data_(size_x, size_y, initial_value)
    {}

    partition_data<T> get_data(direction type) const
    {
            return partition_data<T>(data_, type);
    }
    HPX_DEFINE_COMPONENT_DIRECT_ACTION(partition_server, get_data);

private:
    partition_data<T> data_;
};

}//namespace server
}//namespace grid

#define HPX_REGISTER_PARTITION_SERVER_DECLARATION(...)                      \
    HPX_REGISTER_PARTITION_SERVER_DECLARATION_(__VA_ARGS__)                             \
/**/
#define HPX_REGISTER_PARTITION_SERVER_DECLARATION_(...)                                 \
    HPX_UTIL_EXPAND_(BOOST_PP_CAT(                                            \
        HPX_REGISTER_PARTITION_SERVER_DECLARATION_, HPX_UTIL_PP_NARG(__VA_ARGS__)       \
    )(__VA_ARGS__))                                                           \
/**/

#define HPX_REGISTER_PARITION_SERVER_DECLARATION_1(type)                               \
    HPX_REGISTER_PARTITION_DECLARATION_2(type, type)                             \
/**/
#define HPX_REGISTER_PARTITION_SERVER_DECLARATION_2(type, name)                         \
    HPX_REGISTER_ACTION_DECLARATION(                                          \
        grid::server::partition_server<type>::get_data_action,              \
        BOOST_PP_CAT(__partition_server_get_data_action_, name));                      \
/**/

#define HPX_REGISTER_PARTITION_SERVER(...)                                  \
    HPX_REGISTER_PARTITION_SERVER_(__VA_ARGS__)                                         \
/**/
#define HPX_REGISTER_PARTITION_SERVER_(...)                                             \
    HPX_UTIL_EXPAND_(BOOST_PP_CAT(                                            \
        HPX_REGISTER_PARTITION_SERVER_, HPX_UTIL_PP_NARG(__VA_ARGS__)                   \
    )(__VA_ARGS__))                                                           \
/**/

#define HPX_REGISTER_PARTITION_SERVER_1(type)                                           \
    HPX_REGISTER_PARTITION_SERVER_2(type, type)                                         \
/**/
#define HPX_REGISTER_PARTITION_SERVER_2(type, name)                                     \
    HPX_REGISTER_ACTION(                                                      \
        grid::server::partition_server<type>::get_data_action,            \
        BOOST_PP_CAT(__partition_server_get_data_action_, name));                      \
    typedef ::hpx::components::simple_component<                              \
        grid::server::partition_server<type>                               \
    > BOOST_PP_CAT(__partition_server_, name);                                        \
    HPX_REGISTER_COMPONENT(BOOST_PP_CAT(__partition_server_, name))                     \
/**/

#endif
