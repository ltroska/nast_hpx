#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>
#include <hpx/lcos/gather.hpp>

#include "io/manager.hpp"
#include "stepper/stepper.hpp"

int hpx_main(int argc, char* argv[])
{
    io::config config = io::manager::read_config_from_file("../input.xml");

    stepper::stepper step;
    step.setup(config).wait();

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}
