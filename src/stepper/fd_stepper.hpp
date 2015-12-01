#ifndef STEPPER_STEPPER_HPP
#define STEPPER_STEPPER_HPP

#include "grid/partition.hpp"
#include "internal/cfd_config.hpp"

namespace stepper {

class fd_stepper
{
    public:
        typedef std::vector<std::vector<grid::partition> > space;

        fd_stepper(cfd_config* config);

        space do_work();

    private:
        cfd_config* config;

        RealType dx;
        RealType dy;

        uint numPartitions;
        uint numLocalities;
};


}//namespace stepper
#endif
