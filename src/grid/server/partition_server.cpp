#include "partition_server.hpp"

typedef grid::server::partition_server partition_component;
typedef hpx::components::component<partition_component> partition_server_type;

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_COMPONENT(partition_server_type, partition_component);

HPX_REGISTER_ACTION(grid::server::partition_server::get_p_data_action, partition_server_get_p_data_action);
HPX_REGISTER_ACTION(grid::server::partition_server::get_uv_data_action, partition_server_get_uv_data_action);
HPX_REGISTER_ACTION(grid::server::partition_server::get_fg_data_action, partition_server_get_fg_data_action);
