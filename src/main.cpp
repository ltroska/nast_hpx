#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include "io/config.hpp"
#include "stepper/stepper.hpp"

int hpx_main(boost::program_options::variables_map& vm)
{
    nast_hpx::stepper::stepper step;
        
    const auto cfg_path = vm["cfg"].as<std::string>();
        
    step.setup(nast_hpx::io::config::read_config_from_file(cfg_path.c_str(), hpx::get_locality_id(), hpx::get_num_localities_sync())).wait();
    
    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    using namespace boost::program_options;

    options_description desc_commandline;
    
    desc_commandline.add_options()
    ("cfg", value<std::string>()->default_value("../settings/driven_cavity.xml"),
         "path to config xml");
    
    std::vector<std::string> cfg;
    cfg.push_back("hpx.run_hpx_main!=1");

    
    return hpx::init(desc_commandline, argc, argv, cfg);
}
