#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/util/unwrapped.hpp>

#include <string>

#include "grid/partition.hpp"
#include "io/config_reader.hpp"
#include "stepper/stepper.hpp"

int hpx_main(int argc, char* argv[])
{
 //   if (hpx::get_locality_id() == 0)
 //   {
        cfd_config* config = io::config_reader::read_config_file("input.xml");

        std::vector<hpx::id_type> localities = hpx::find_all_localities();
        uint num_localities = localities.size();

        std::vector<stepper::stepper> steppers;
        std::vector<hpx::future<uint>> futures;

      //  for (uint i = 0; i < num_localities; i++)
      //  {
            //hpx::cout << "stepper on " << i << hpx::endl << hpx::flush;
            stepper::stepper step = stepper::stepper();
            step.setup(32, 32, 1, 1, 2, 2);
            hpx::this_thread::sleep_for(boost::chrono::seconds(4));
            step.do_work(3,3,3,3,3,3).wait();
            //futures.push_back(steppers[i].setup(16, 16, 1, 1, 2 , 2));
      //  }

      //  hpx::when_each(hpx::util::unwrapped([&](uint t) {
       //     steppers[t].do_work(3,3,3,3,3,3);
      //  }),
      //  futures);

       // hpx::cout << "started do work ... " << hpx::endl << hpx::flush;

   // }
    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    // Initialize and run HPX, this executes hpx_main above.
    return hpx::init(argc, argv);
}
