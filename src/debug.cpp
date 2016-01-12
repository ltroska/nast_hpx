#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/util/unwrapped.hpp>

#include <string>
#include <algorithm>

#include "grid/partition.hpp"
#include "io/config_reader.hpp"
#include "stepper/stepper.hpp"

RealType compute_delta_t(RealType u_max, RealType v_max, RealType dx, RealType dy, RealType Re, RealType tau) {
    return tau*std::min(Re/2. * 1./( 1./(dx*dx) + 1./(dy*dy)), std::min(dx/abs(u_max), dy/abs(v_max)) );
}

int hpx_main(int argc, char* argv[])
{
        cfd_config config = io::config_reader::read_config_file("input.xml");

        stepper::stepper step = stepper::stepper();
        step.setup(config);

        if(hpx::get_locality_id() == 0) {
            std::vector<hpx::id_type> localities = hpx::find_all_localities();
            uint num_localities = localities.size();

            std::vector<stepper::stepper> steppers;
            std::vector<hpx::future<uint>> futures;

            for(int i = 0; i < num_localities; i++)
            {
                steppers.push_back(stepper::stepper(hpx::find_from_basename(stepper::server::stepper_basename, i)));
            }

            RealType t = 0;
            RealType delta_t = config.deltaT;
            RealType dx = config.xLength/config.iMax;
            RealType dy = config.yLength/config.jMax;

            grid::vector_cell max_uv = grid::vector_cell(0, 0);

            for(uint n = 0; t < config.tEnd; t += delta_t, n++) {
                delta_t = compute_delta_t(max_uv.u, max_uv.v, dx, dy, config.Re, config.tau);

                for(stepper::stepper stepper : steppers)
                {
                    futures.push_back(stepper.update_delta_t(delta_t));
                }

                hpx::wait_all(futures);
                futures.clear();

                for(stepper::stepper stepper : steppers)
                {
                    futures.push_back(stepper.set_velocity_on_boundary());
                }

                hpx::wait_all(futures);
                futures.clear();

                for(auto stepper : steppers)
                {
                    futures.push_back(stepper.compute_fg());
                }

                hpx::wait_all(futures);

                futures.clear();

                for(auto stepper : steppers)
                {
                    futures.push_back(stepper.set_rhs());
                }

                hpx::wait_all(futures);
                futures.clear();


                //----------------- SOR LOOP -----------------//
                double residual;
                uint i = 0;

                do {
                    i++;

                    for(auto stepper : steppers)
                    {
                        futures.push_back(stepper.set_pressure_on_boundary());
                    }

                    hpx::wait_all(futures);
                    futures.clear();


                    for(auto stepper : steppers)
                    {
                        futures.push_back(stepper.sor_cycle());
                    }

                    hpx::wait_all(futures);
                    futures.clear();

                    hpx::future<RealType> res = hpx::make_ready_future(0.0);
                    for(auto stepper : steppers)
                            res = hpx::lcos::local::dataflow(
                                    [](hpx::future<RealType> next_summand, hpx::future<RealType> prev_sum)
                                         {
                                            return prev_sum.get() + next_summand.get();
                                         }
                                    , stepper.get_residual()
                                    , res
                                    );
                    residual = res.get();
                } while (residual > config.eps*config.eps);


                hpx::future<grid::vector_cell> new_max_uv = hpx::make_ready_future(max_uv);
                for(auto stepper : steppers)
                {
                    new_max_uv = hpx::lcos::local::dataflow(
                        [](hpx::future<grid::vector_cell> max_uv, hpx::future<grid::vector_cell> curr_max_uv)
                            {
                                grid::vector_cell loc_max_uv = max_uv.get();
                                grid::vector_cell curr_uv = curr_max_uv.get();
                                return grid::vector_cell( (loc_max_uv.u >= curr_uv.u) ? loc_max_uv.u : curr_uv.u,
                                                            (loc_max_uv.v >= curr_uv.v) ? loc_max_uv.v : curr_uv.v);
                            }
                            , stepper.update_velocities()
                            , new_max_uv
                            );
                }

                max_uv = new_max_uv.get();

                if (max_uv.u < 2)
                    max_uv.u = 2;

                if (max_uv.v < 2)
                    max_uv.v = 2;


                hpx::cout << "step " << n << " dt " << delta_t << " residual " << residual << " iterations " << i << " max uv " << max_uv.u << " " << max_uv.v << hpx::endl << hpx::flush;

                for(auto stepper : steppers)
                {
                    futures.push_back(stepper.do_work(3,3,3,3,3,3));
                }

                hpx::wait_all(futures);
                futures.clear();
            }

        }

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    // Initialize and run HPX, this executes hpx_main above.
    return hpx::init(argc, argv);
}
