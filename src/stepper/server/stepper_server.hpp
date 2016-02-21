#ifndef STEPPER_SERVER_STEPPER_HPP
#define STEPPER_SERVER_STEPPER_HPP

#include <hpx/include/components.hpp>

#include "io/config.hpp"
#include "computation/parameters.hpp"
#include "computation/strategy.hpp"
#include "grid/partition.hpp"

namespace stepper { namespace server {

char const* stepper_basename = "/cfd_hpx/stepper/";

struct HPX_COMPONENT_EXPORT stepper_server
    : hpx::components::component_base<stepper_server>
{
    public:
        stepper_server() {}
        stepper_server(uint num_localities);
        stepper_server(uint num_localities, io::config const& cfg);

    protected:
        template<typename T>
        void print_grid(std::vector<grid::partition<T> > const& grid, const std::string message = "") const;

        uint get_index(uint k, uint l) const;

    private:
        void initialize_parameters();
        void initialize_grids();

        io::config c;
        computation::parameters params;
        computation::strategy* strategy;

        uint num_localities, num_localities_x, num_localities_y;
        uint num_partitions_x, num_partitions_y;

        index_grid_type index_grid;
        vector_grid_type uv_grid;



};

}//namespace server
}//namespace stepper
#endif
