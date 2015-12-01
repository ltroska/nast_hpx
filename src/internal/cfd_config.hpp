#ifndef INTERNAL_CFD_CONFIG_HPP
#define INTERNAL_CFD_CONFIG_HPP

#include "types.hpp"

 struct cfd_config {
        uint iMax;
        uint jMax;
        RealType xLength;
        RealType yLength;
        RealType Re;
        RealType omega;
        RealType tau;
        RealType eps;
        RealType alpha;
        uint iterMax;
};

#endif
