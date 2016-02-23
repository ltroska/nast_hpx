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

        void do_timestep(RealType dt);

    protected:
        template<typename T>
        void print_grid(std::vector<grid::partition<T> > const& grid, const std::string message = "") const;

        uint get_index(uint k, uint l) const;

    private:
        void initialize_parameters();
        void initialize_grids();

        RealType compute_new_dt(RealType u_max, RealType v_max);

        io::config c;
        computation::parameters params;
        computation::strategy* strategy;

        uint num_localities, num_localities_x, num_localities_y;

        index_grid_type index_grid;
        vector_grid_type uv_grid, fg_grid;
        scalar_grid_type p_grid, rhs_grid;

};

}//namespace server
}//namespace stepper
#endif
