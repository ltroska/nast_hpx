#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/runtime/actions/plain_action.hpp>

#include<string>

#include "grid/partition.hpp"
#include "io/config_reader.hpp"

void test(grid::partition const& left, grid::partition const& right) {
    hpx::cout << left.get_data(grid::left_partition).get();
}

HPX_PLAIN_ACTION(test, test_action);

int hpx_main(boost::program_options::variables_map& vm)
{
    std::string config_path = vm["config"].as<std::string>();
    cfd_config* config = io::config_reader::read_config_file(config_path.c_str());

    std::vector<hpx::id_type> localities = hpx::find_all_localities();
    uint nl = localities.size();

    hpx::cout << "#localities = " << nl << "\n";

    grid::partition par(hpx::find_here(), 3,3, 1);
    grid::partition par2(hpx::find_here(), 3,3, 10);


    test_action act;
    hpx::async(act, hpx::find_here(), par, par2);

    // Initiate shutdown of the runtime system.
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
