#include "config.hpp"

#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>

#include "pugixml/pugixml.hpp"

namespace nast_hpx { namespace io {
    /// Methods reads the simulation configuration from the given file
    /// and returns a corresponding config object.
    config config::read_config_from_file(const char *path, std::size_t rank, std::size_t num_localities)
    {
        config cfg;
        cfg.num_localities = num_localities;
        cfg.rank = rank;
        
        if (cfg.num_localities == 2)
        {
            cfg.num_localities_x = 2;
            cfg.num_localities_y = 1;
        }
        else
        {
            cfg.num_localities_x = static_cast<uint> (sqrt(cfg.num_localities));
            cfg.num_localities_y = cfg.num_localities_x;
        }
    
        pugi::xml_document doc;

        pugi::xml_parse_result result = doc.load_file(path);

        if (!result) {
            std::cerr << "Error loading file: " << path << std::endl;
            std::exit(1);
        }

        pugi::xml_node config_node = doc.child("SimulationConfig");
        if (config_node == NULL) {
            std::cerr 
                << "Error: A simulation configuration must be defined!" 
                << std::endl;
            std::exit(1);
        }

        if(config_node.child("iMax") != NULL) {
            cfg.i_max = config_node.child("iMax").first_attribute().as_int();
        }
        else {
            std::cerr << "Error: iMax not set!" << std::endl;
            std::exit(1);
        }

        if(config_node.child("jMax") != NULL)
        {
            cfg.j_max = config_node.child("jMax").first_attribute().as_int();
        }
        else {
            std::cerr << "Error: jMax not set!" << std::endl;
            std::exit(1);
        }

        if(config_node.child("iRes") != NULL)
        {
            cfg.i_res = config_node.child("iRes").first_attribute().as_int();
        }
        else {
            std::cerr << "Error: iRes not set!" << std::endl;
            std::exit(1);
        }

        if(config_node.child("jRes") != NULL)
        {
            cfg.j_res = config_node.child("jRes").first_attribute().as_int();
        }
        else {
            std::cerr << "Error: jRes not set!" << std::endl;
            std::exit(1);
        }
        
        cfg.num_local_partitions_x = cfg.i_res;
        cfg.num_local_partitions_y = cfg.j_res;
        cfg.num_local_partitions = cfg.num_local_partitions_x * cfg.num_local_partitions_y;
        
        cfg.cells_x_per_partition = (cfg.i_max + 2) / cfg.num_localities_x / cfg.num_local_partitions_x;
        
        if (cfg.cells_x_per_partition * cfg.num_localities_x * cfg.num_local_partitions_x != cfg.i_max + 2)
        {
            std::cerr << "Error: localities_x * num_local_partitions_x does not divide i_max + 2 evenly!" << std::endl;
            std::exit(1);
        }        
        
        cfg.cells_y_per_partition = (cfg.j_max + 2) / cfg.num_localities_y / cfg.num_local_partitions_y;
        
        if (cfg.cells_y_per_partition * cfg.num_localities_y * cfg.num_local_partitions_y != cfg.j_max + 2)
        {
            std::cerr << "Error: localities_y * num_local_partitions_y does not divide j_max + 2 evenly!" << std::endl;
            std::exit(1);
        }

        if(config_node.child("xLength") != NULL)
        {
            cfg.x_length =
                config_node.child("xLength").first_attribute().as_double();
        }
        else {
            std::cerr << "Error: xLength not set!" << std::endl;
            std::exit(1);
        }

        if(config_node.child("yLength") != NULL)
        {
            cfg.y_length =
                config_node.child("yLength").first_attribute().as_double();
        }
        else {
            std::cerr << "Error: yLength not set!" << std::endl;
            std::exit(1);
        }

        if(config_node.child("Re") != NULL)
        {
            cfg.re = config_node.child("Re").first_attribute().as_double();
        }
        else {
            std::cerr << "Error: Re not set!" << std::endl;
            std::exit(1);
        }
           

        if(config_node.child("Pr") != NULL)
        {
            cfg.pr = config_node.child("Pr").first_attribute().as_double();
        }
        else
        {
            cfg.pr = 0;
        }

        if(config_node.child("omega") != NULL)
        {
            cfg.omega =
                config_node.child("omega").first_attribute().as_double();
        }
        else
        {
            std::cerr << "Error: Omega not set!" << std::endl;
            std::exit(1);
        }
        
        if(config_node.child("tau") != NULL)
        {
            cfg.tau = config_node.child("tau").first_attribute().as_double();
        }
        else
        {
            cfg.tau = 0.5;
        }

        if(config_node.child("eps") != NULL)
        {
            cfg.eps = config_node.child("eps").first_attribute().as_double();
            cfg.eps_sq = cfg.eps * cfg.eps;
        }
        else
        {
            std::cerr << "Error: Eps not set!" << std::endl;
            std::exit(1);
        }
        
        if(config_node.child("alpha") != NULL)
        {
            cfg.alpha =
                config_node.child("alpha").first_attribute().as_double();
        }
        else
        {
            cfg.alpha = 0.9;
        }

        if(config_node.child("beta") != NULL)
        {
            cfg.beta = config_node.child("beta").first_attribute().as_double();
        }
        else
        {
            cfg.beta = 0;
        }

        if(config_node.child("iterMax") != NULL)
        {
            cfg.iter_max =
                config_node.child("iterMax").first_attribute().as_int();
        }
        else
        {
            std::cerr << "Error: iterMax not set!" << std::endl;
            std::exit(1);
        }

        if(config_node.child("tEnd") != NULL)
        {
            cfg.t_end = config_node.child("tEnd").first_attribute().as_double();
        }
        else
        {
            std::cerr << "Error: tEnd not set!" << std::endl;
            std::exit(1);
        }

        if(config_node.child("dt") != NULL)
        {
            cfg.dt = config_node.child("dt").first_attribute().as_double();
        }
        else
            cfg.dt = 0.01;

        if(config_node.child("outputSkipSize") != NULL)
        {
            cfg.output_skip_size =
                config_node.child("outputSkipSize").first_attribute().as_int();
        }
        else
        {
            cfg.output_skip_size = 1;
        }

        if(config_node.child("subIterations") != NULL)
        {
            cfg.sub_iterations =
                config_node.child("subIterations").first_attribute().as_int();
        }
        else
        {
            cfg.sub_iterations = 1;
        }

        if(config_node.child("withForEach") != NULL)
        {
            cfg.wfe =
                config_node.child("withForEach").first_attribute().as_int();
        }
        else
        {
            cfg.wfe = 1;
        }

        if(config_node.child("vtk") != NULL)
        {
            cfg.vtk =
                (config_node.child("vtk").first_attribute().as_int() == 1);
        }
        else
        {
            cfg.vtk = 0;
        }

        if(config_node.child("GX") != NULL)
        {
            cfg.gx = config_node.child("GX").first_attribute().as_double();
        }
        else
        {
            cfg.gx = 0;
        }

        if(config_node.child("GY") != NULL)
        {
            cfg.gy = config_node.child("GY").first_attribute().as_double();
        }
        else
        {
            cfg.gy = 0;
        }

        if(config_node.child("UO") != NULL)
        {
            cfg.u_bnd.top =
                config_node.child("UO").first_attribute().as_double();
        }
        else
        {
            cfg.u_bnd.top = 0;
        }

        if(config_node.child("UB") != NULL)
        {
            cfg.u_bnd.bottom =
                config_node.child("UB").first_attribute().as_double();
        }
        else
        {
            cfg.u_bnd.bottom = 0;
        }

        if(config_node.child("UL") != NULL)
        {
            cfg.u_bnd.left =
                config_node.child("UL").first_attribute().as_double();
        }
        else
        {
            cfg.u_bnd.left = 0;
        }

        if(config_node.child("UR") != NULL)
        {
            cfg.u_bnd.right =
                config_node.child("UR").first_attribute().as_double();
        }
        else
        {
            cfg.u_bnd.right = 0;
        }

        if(config_node.child("VO") != NULL)
        {
            cfg.v_bnd.top =
                config_node.child("VO").first_attribute().as_double();
        }
        else
        {
            cfg.v_bnd.top = 0;
        }

        if(config_node.child("VB") != NULL)
        {
            cfg.v_bnd.bottom =
                config_node.child("VB").first_attribute().as_double();
        }
        else
        {
            cfg.v_bnd.bottom = 0;
        }

        if(config_node.child("VL") != NULL)
        {
            cfg.v_bnd.left =
                config_node.child("VL").first_attribute().as_double();
        }
        else
        {
            cfg.v_bnd.left = 0;
        }

        if(config_node.child("VR") != NULL)
        {
            cfg.v_bnd.right =
                config_node.child("VR").first_attribute().as_double();
        }
        else
        {
            cfg.v_bnd.right = 0;
        }

        if(config_node.child("TO") != NULL)
        {
            cfg.temp_bnd.top =
                config_node.child("TO").first_attribute().as_double();
        }
        else
        {
            cfg.temp_bnd.top = 0;
        }

        if(config_node.child("TB") != NULL)
        {
            cfg.temp_bnd.bottom =
                config_node.child("TB").first_attribute().as_double();
        }
        else
        {
            cfg.temp_bnd.bottom = 0;
        }

        if(config_node.child("TL") != NULL)
        {
            cfg.temp_bnd.left =
                config_node.child("TL").first_attribute().as_double();
        }
        else
        {
            cfg.temp_bnd.left = 0;
        }

        if(config_node.child("TR") != NULL)
        {
            cfg.temp_bnd.right =
                config_node.child("TR").first_attribute().as_double();
        }
        else
        {
            cfg.temp_bnd.right = 0;
        }

        if(config_node.child("TI") != NULL)
        {
            cfg.ti = config_node.child("TI").first_attribute().as_double();
        }
        else
        {
            cfg.ti = 0;
        }

        if(config_node.child("WTL") != NULL)
        {
            cfg.temp_data_type.left =
                config_node.child("WTL").first_attribute().as_double();
        }
        else
        {
            cfg.temp_data_type.left = -1;
        }

        if(config_node.child("WTR") != NULL)
        {
            cfg.temp_data_type.right =
                config_node.child("WTR").first_attribute().as_double();
        }
        else
        {
            cfg.temp_data_type.right = -1;
        }

        if(config_node.child("WTB") != NULL)
        {
            cfg.temp_data_type.bottom =
                config_node.child("WTB").first_attribute().as_double();
        }
        else
        {
            cfg.temp_data_type.bottom = -1;
        }

        if(config_node.child("WTT") != NULL)
        {
            cfg.temp_data_type.top =
                config_node.child("WTT").first_attribute().as_double();
        }
        else
        {
            cfg.temp_data_type.top = -1;
        }

        if(config_node.child("WL") != NULL)
        {
            cfg.data_type.left =
                config_node.child("WL").first_attribute().as_double();
        }
        else
        {
            cfg.data_type.left = 1;
        }

        if(config_node.child("WR") != NULL)
        {
            cfg.data_type.right =
                config_node.child("WR").first_attribute().as_double();
        }
        else
        {
            cfg.data_type.right = 1;
        }

        if(config_node.child("WB") != NULL)
        {
            cfg.data_type.bottom =
                config_node.child("WB").first_attribute().as_double();
        }
        else
        {
            cfg.data_type.bottom = 1;
        }

        if(config_node.child("WT") != NULL)
        {
            cfg.data_type.top =
                config_node.child("WT").first_attribute().as_double();
        }
        else
        {
            cfg.data_type.top = 1;
        }

        if(config_node.child("deltaVec") != NULL)
        {
            cfg.delta_vec =
                config_node.child("deltaVec").first_attribute().as_double();
        }
        else
        {
            cfg.delta_vec = 0;
        }

        if(config_node.child("gridFile") != NULL)
        {
            cfg.with_flag_grid = true;

            std::ifstream file(
                config_node.child("gridFile").first_attribute().as_string());
            
            cfg.num_fluid_cells = 0;

            std::size_t flag_res_x = cfg.num_local_partitions_x * cfg.cells_x_per_partition + 2;
            std::size_t flag_res_y = cfg.num_local_partitions_y * cfg.cells_y_per_partition + 2;
            
          //  std::cout << "resx " << flag_res_x << " resy " << flag_res_y << std::endl;

            cfg.flag_grid.resize(flag_res_x * flag_res_y);
                                    
          //  std::cout << "s " << cfg.flag_grid.size() << std::endl;
            std::size_t i = 0;
            std::size_t j = cfg.j_max + 1;
            
            std::size_t start_i = (rank % cfg.num_localities_x) * cfg.num_local_partitions_x * cfg.cells_x_per_partition;
            std::size_t end_i = start_i + cfg.num_local_partitions_x * cfg.cells_x_per_partition;            
            
            std::size_t start_j = (rank / cfg.num_localities_x) * cfg.num_local_partitions_y * cfg.cells_y_per_partition;
            std::size_t end_j = start_j + cfg.num_local_partitions_y * cfg.cells_y_per_partition;
            
           // std::cout <<  " starti " << start_i << " endi " << end_i << std::endl;
           // std::cout <<  " startj " << start_j << " endj " << end_j << std::endl;
            
            std::size_t offset_x = 0;
            std::size_t offset_y = flag_res_y - 1;
            
            if (start_i == 0)
            {
                for (std::size_t j = 0; j < flag_res_y; ++j)
                    cfg.flag_grid[j * flag_res_x + 0] = std::bitset<6>("000000");

                ++offset_x;
            }
            else
                --start_i;
            
            if (start_j == 0)
            {
                for (std::size_t i = 0; i < flag_res_x; ++i)
                    cfg.flag_grid[0 * flag_res_x + i] = std::bitset<6>("000000");

            }
            else
                --start_j; 

                
            if (end_i == cfg.i_max + 2)
            {
                for (std::size_t j = 0; j < flag_res_y; ++j)
                    cfg.flag_grid[j * flag_res_x + flag_res_x - 1] = std::bitset<6>("000000");
                    
                --end_i;
            }  
              
            if (end_j == cfg.j_max + 2)
            {
                for (std::size_t i = 0; i < flag_res_x; ++i)
                    cfg.flag_grid[(flag_res_y - 1) * flag_res_x + i] = std::bitset<6>("000000");
                    
                --end_j;
                --offset_y;
            }
                
          //  std::cout <<  " starti " << start_i << " endi " << end_i << std::endl;
            //std::cout <<  " startj " << start_j << " endj " << end_j << std::endl;

                
            std::size_t insert_i = 0;
            std::size_t insert_j = offset_y;

            while (true)
            {
                std::string line;
                std::getline(file, line);

                if (!file.good())
                    break;

                std::stringstream iss(line);
                i = 0;
                insert_i = offset_x;
                while (true)
                {
                    std::string cell_val;
                    std::getline(iss, cell_val, ',');
                    
                    std::bitset<6> flag(std::stoi(cell_val));

                    if (flag.test(4))
                        cfg.num_fluid_cells++;
                    
                    if (i >= start_i && i <= end_i && j >= start_j && j <= end_j)
                    {
                        cfg.flag_grid[insert_j * flag_res_x + insert_i] = flag;
                        ++insert_i;
                    }
                    
                    if (!iss.good())
                        break;
                    ++i;
                }
                
                if (j >= start_j && j <= end_j)
                    --insert_j;

                --j;
            }
            
         /*   for (uint j = flag_res_y - 1; j < flag_res_y; --j)
            {
                for (uint i = 0; i < flag_res_x; ++i)
                    std::cout << cfg.flag_grid[j * flag_res_x + i].to_ulong() << " ";
                std::cout << "\n";
                    
            }*/
            
            std::cout << std::endl;
        }
        else
        {
            std::cerr << "Error: geometry file not set!" << std::endl;
            std::exit(1);
        }

        if(config_node.child("initialUVFile") != NULL)
        {
            cfg.with_initial_uv_grid = true;

            std::ifstream file(
                config_node.child("initialUVFile")
                    .first_attribute().as_string());

            while (true)
            {
                std::string line;
                std::getline(file, line);

                if (!file.good())
                    break;

                std::stringstream iss(line);

                while (true)
                {
                    std::string cell_val;
                    std::getline(iss, cell_val, ',');

                    std::stringstream issc(cell_val);

                    std::string u;
                    std::string v;
                    std::getline(issc, u, '/');
                    std::getline(issc, v, '/');

                    cfg.initial_uv_grid.push_back(
                        std::pair<Real, Real>(
                            std::stod(u), std::stod(v)));

                    if (!iss.good())
                        break;
                }
            }


        }
        else
        {
            cfg.with_initial_uv_grid = false;
        }

        return cfg;
    }
        
}
}