#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/runtime/actions/plain_action.hpp>

#include<string>

#include "grid/partition.hpp"
#include "io/config_reader.hpp"
#include "stepper/stepper.hpp"

int hpx_main(boost::program_options::variables_map& vm)
{
    std::string config_path = vm["config"].as<std::string>();
    cfd_config* config = io::config_reader::read_config_file(config_path.c_str());

    stepper::stepper stepper = stepper::stepper();

   // stepper.do_work();

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    namespace po = boost::program_options;

    po::options_description desc_commandline;
    desc_commandline.add_options()
        ("help,h", "print help message")
        ("results,r", "print generated results (default: false)")
        ("config,c", po::value<std::string>(),
         "path to config file")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc_commandline), vm);
    po::notify(vm);

    if(vm.count("help"))
    {
        std::cout << desc_commandline << std::endl;
        return EXIT_FAILURE;
    }

    if(!vm.count("config"))
    {
        std::cerr << "The --config/--c option is required." << std::endl;
        return EXIT_FAILURE;
    }
    // Initialize and run HPX, this executes hpx_main above.
    return hpx::init(desc_commandline, argc, argv);
}
