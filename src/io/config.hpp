#ifndef IO_CONFIG_HPP
#define IO_CONFIG_HPP

#include "util/types.hpp"

namespace io {

struct config
{
        uint i_max;
        uint j_max;
        uint i_res;
        uint j_res;
        RealType x_length;
        RealType y_length;
        RealType re;
        RealType omega;
        RealType tau;
        RealType eps;
        RealType eps_sq;
        RealType alpha;
        RealType t_end;
        RealType delta_t;
        uint output_skip_size;
        uint sub_iterations;
        uint iter_max;

        uint wfe;

        template <typename Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            ar & i_max & j_max & i_res & j_res & x_length & y_length & re & omega & tau & eps & alpha & t_end & delta_t & output_skip_size & iter_max;
        }

        friend std::ostream& operator<<(std::ostream& os, config const& config)
        {
            os << "config: {iMax = " << config.i_max << ", jMax = "<< config.j_max
                << ", iRes = " << config.i_res << ", jRes = " << config.j_res << ", xLength = " << config.x_length
                << ", yLength = " << config.y_length << ", Re = " << config.re << ", omega = " << config.omega
                << ", tau = " << config.tau << " eps = " << config.eps << ", alpha = " << config.alpha
                << ", outputSkipSize = " << config.output_skip_size << ", iterMax = " << config.iter_max << ", wfe = " << config.wfe << "}";
            return os;
        }
};

} //namespace

#endif

