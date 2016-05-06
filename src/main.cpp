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

    auto a1 = scalar_cell(1);
    auto a2 = 1.;
    
    a1 = 3;
    
    auto b1 = vector_cell(1., 1.);
    auto b2 = std::make_pair(1., 1.);
    
    std::cout << sizeof(a1) << " " << sizeof(a2) << std::endl;
    std::cout << sizeof(b1) << " " << sizeof(b2) << std::endl;
    std::cout << sizeof(grid::partition_data<std::bitset<5> >) << " " << sizeof(std::vector<std::bitset<5> >) << std::endl;

    stepper::stepper step;
        
    step.setup(io::config::read_config_from_file(argv[1])).wait();

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}
