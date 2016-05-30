#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include "io/config.hpp"
#include "stepper/stepper.hpp"

#include <hpx/lcos/when_each.hpp>

int hpx_main(boost::program_options::variables_map& vm)
{
   nast_hpx::stepper::stepper step;
    
  /*  std::vector<hpx::future<int> > a;
    a.push_back(hpx::make_ready_future(1)); 
    a.push_back(hpx::make_ready_future(2)); 
    
    auto f1 = [](hpx::future<int> f) {std::cout << "got " << f.get() << " on first " << std::endl;};
    auto f2 = [](int idx, hpx::future<int> f) {std::cout << "got " << f.get() << " with index " << idx << "on second" << std::endl;};
    
    hpx::when_each_n(f2, a.begin(), 1);*/
        
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