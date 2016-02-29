#include "partition_server.hpp"
#include "util/cell.hpp"

HPX_REGISTER_COMPONENT_MODULE();

/*
*@TODO: Why does this work???
*/

//HPX_REGISTER_PARTITION_SERVER_DECLARATION(scalar_cell);
//HPX_REGISTER_PARTITION_SERVER_DECLARATION(vector_cell);
HPX_REGISTER_PARTITION_SERVER(scalar_cell);
HPX_REGISTER_PARTITION_SERVER(vector_cell);
