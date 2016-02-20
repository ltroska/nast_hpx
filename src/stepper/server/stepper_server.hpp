#ifndef STEPPER_SERVER_STEPPER_HPP
#define STEPPER_SERVER_STEPPER_HPP

#include <hpx/include/components.hpp>

#include "io/config.hpp"
#include "solver/parameters.hpp"
#include "solver/solver.hpp"
#include "grid/partition.hpp"

namespace stepper { namespace server {

char const* stepper_basename = "/cfd_hpx/stepper/";

struct HPX_COMPONENT_EXPORT stepper_server
    : hpx::components::component_base<stepper_server>
{
    protected:
        typedef std::vector<grid::partition> grid_type;
        typedef std::vector<std::pair<uint, uint> > index_grid_type;

    public:
        stepper_server() {}
        stepper_server(uint num_localities);
        stepper_server(uint num_localities, io::config const& cfg);

    protected:
        void print_grid(grid_type const& grid, const std::string message = "") const;

        uint get_index(uint k, uint l) const;

    private:
        io::config c;
        solver::parameters params;
        solver::solver* solv;

        uint num_localities, num_localities_x, num_localities_y;
        uint num_partitions_x, num_partitions_y;

        index_grid_type index_grid;
        grid_type u_grid, v_grid;



};

}//namespace server
}//namespace stepper
#endif
