#ifndef INTERNAL_CFD_CONFIG_HPP
#define INTERNAL_CFD_CONFIG_HPP

#include "types.hpp"

struct cfd_config
{
        uint iMax;
        uint jMax;
        uint iRes;
        uint jRes;
        RealType xLength;
        RealType yLength;
        RealType Re;
        RealType omega;
        RealType tau;
        RealType eps;
        RealType alpha;
        RealType tEnd;
        RealType deltaT;
        uint iterMax;

        friend std::ostream& operator<<(std::ostream& os, cfd_config const& config)
        {
            os << "config: {iMax = " << config.iMax << ", jMax = "<< config.jMax
                << ", iRes = " << config.iRes << ", jRes = " << config.jRes << ", xLength = " << config.xLength
                << ", yLength = " << config.yLength << ", Re = " << config.Re << ", omega = " << config.omega
                << ", tau = " << config.tau << " eps = " << config.eps << ", alpha = " << config.alpha
                << ", iterMax = " << config.iterMax << "}" << std::endl;
            return os;
        }
};

#endif
