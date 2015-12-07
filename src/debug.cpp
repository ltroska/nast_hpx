#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/runtime/actions/plain_action.hpp>

#include <string>

#include "grid/partition.hpp"
#include "io/config_reader.hpp"
#include "stepper/stepper.hpp"

int hpx_main(int argc, char* argv[])
{
    cfd_config* config = io::config_reader::read_config_file("input.xml");

    std::vector<hpx::id_type> localities = hpx::find_all_localities();
    uint num_localities = localities.size();

    stepper::stepper step = stepper::stepper();
    hpx::future<uint> result = step.do_work(2, 2, 2, 2, 0.1, 0.1);

    uint i = result.get();

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    // Initialize and run HPX, this executes hpx_main above.
    return hpx::init(argc, argv);
}
