#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include "io/config.hpp"
#include "stepper/stepper.hpp"

int hpx_main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: ./main <input.xml>" << std::endl;
        exit(0);
    }

    stepper::stepper step;
        
    step.setup(io::config::read_config_from_file(argv[1])).wait();

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}
