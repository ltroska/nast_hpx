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

        cfg.num_localities_x = static_cast<uint> (std::cbrt(cfg.num_localities));
        cfg.num_localities_y = cfg.num_localities_x;
        cfg.num_localities_z = cfg.num_localities_x;

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

        if(config_node.child("kMax") != NULL)
        {
            cfg.k_max = config_node.child("kMax").first_attribute().as_int();
        }
        else {
            std::cerr << "Error: kMax not set!" << std::endl;
            std::exit(1);
        }

        if(config_node.child("iRes") != NULL)
        {
            cfg.num_x_blocks = config_node.child("iRes").first_attribute().as_int();
        }
        else {
            std::cerr << "Error: iRes not set!" << std::endl;
            std::exit(1);
        }

        if(config_node.child("jRes") != NULL)
        {
            cfg.num_y_blocks = config_node.child("jRes").first_attribute().as_int();
        }
        else {
            std::cerr << "Error: jRes not set!" << std::endl;
            std::exit(1);
        }

        if(config_node.child("kRes") != NULL)
        {
            cfg.num_z_blocks = config_node.child("kRes").first_attribute().as_int();
        }
        else {
            std::cerr << "Error: kRes not set!" << std::endl;
            std::exit(1);
        }

        cfg.cells_x_per_block = (cfg.i_max + 2) / cfg.num_localities_x / cfg.num_x_blocks;

        if (cfg.cells_x_per_block * cfg.num_localities_x * cfg.num_x_blocks != cfg.i_max + 2)
        {
            std::cerr << "Error: localities_x * num_x_blocks does not divide i_max + 2 evenly!" << std::endl;
            std::cerr << "localities_x = " << cfg.num_localities_x << ", num_x_blocks = " << cfg.num_x_blocks << ", i_max + 2 = " << cfg.i_max + 2 << std::endl;
            std::exit(1);
        }

        cfg.cells_y_per_block = (cfg.j_max + 2) / cfg.num_localities_y / cfg.num_y_blocks;

        if (cfg.cells_y_per_block * cfg.num_localities_y * cfg.num_y_blocks != cfg.j_max + 2)
        {
            std::cerr << "Error: localities_y * num_y_blocks does not divide j_max + 2 evenly!" << std::endl;
            std::cerr << "localities_y = " << cfg.num_localities_y << ", num_y_blocks = " << cfg.num_y_blocks << ", j_max + 2 = " << cfg.j_max + 2 << std::endl;
            std::exit(1);
        }

        cfg.cells_z_per_block = (cfg.k_max + 2) / cfg.num_localities_z / cfg.num_z_blocks;

        if (cfg.cells_z_per_block * cfg.num_localities_z * cfg.num_z_blocks != cfg.k_max + 2)
        {
            std::cerr << "Error: localities_z * num_z_blocks does not divide k_max + 2 evenly!" << std::endl;
            std::cerr << "localities_z = " << cfg.num_localities_z << ", num_z_blocks = " << cfg.num_z_blocks << ", k_max + 2 = " << cfg.k_max + 2 << std::endl;
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

        if(config_node.child("zLength") != NULL)
        {
            cfg.z_length =
                config_node.child("zLength").first_attribute().as_double();
        }
        else {
            std::cerr << "Error: zLength not set!" << std::endl;
            std::exit(1);
        }

        cfg.dx = cfg.x_length / cfg.i_max;
        cfg.dy = cfg.y_length / cfg.j_max;
        cfg.dz = cfg.z_length / cfg.k_max;

        cfg.dx_sq = std::pow(cfg.dx, 2);
        cfg.dy_sq = std::pow(cfg.dy, 2);
        cfg.dz_sq = std::pow(cfg.dz, 2);

        cfg.over_dx = 1. / cfg.dx;
        cfg.over_dy = 1. / cfg.dy;
        cfg.over_dz = 1. / cfg.dz;

        cfg.over_dx_sq = 1. / cfg.dx_sq;
        cfg.over_dy_sq = 1. / cfg.dy_sq;
        cfg.over_dz_sq = 1. / cfg.dz_sq;

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

        cfg.part1 = 1. - cfg.omega;
        cfg.part2 = cfg.omega * cfg.dx_sq * cfg.dy_sq * cfg.dz_sq / (2. * (cfg.dx_sq * cfg.dy_sq + cfg.dx_sq * cfg.dz_sq + cfg.dy_sq * cfg.dz_sq));
        cfg.factor_jacobi = cfg.dx_sq * cfg.dy_sq / (2. * (cfg.dx_sq + cfg.dy_sq));

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

        if(config_node.child("maxTimesteps") != NULL)
        {
            cfg.max_timesteps = config_node.child("maxTimesteps").first_attribute().as_int();
        }
        else
        {
            cfg.max_timesteps = 0;
        }

        if(config_node.child("dt") != NULL)
        {
            cfg.initial_dt = config_node.child("dt").first_attribute().as_double();
        }
        else
            cfg.initial_dt = 0.01;

        if(config_node.child("subIterations") != NULL)
        {
            cfg.sub_iterations =
                config_node.child("subIterations").first_attribute().as_int();
        }
        else
        {
            cfg.sub_iterations = 1;
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

        if(config_node.child("GZ") != NULL)
        {
            cfg.gz = config_node.child("GZ").first_attribute().as_double();
        }
        else
        {
            cfg.gz = 0;
        }

        if(config_node.child("BoundaryConditions") != NULL)
        {
            auto bc_node = config_node.child("BoundaryConditions");

            if (bc_node.child("Top") != NULL)
            {
                auto node = bc_node.child("Top");
                std::string type = node.attribute("type").value();

                if (type == "noslip")
                    cfg.bnd_condition.top_type = noslip;
                else if (type == "slip")
                    cfg.bnd_condition.top_type = slip;
                else if (type == "instream")
                    cfg.bnd_condition.top_type = instream;
                else if (type == "outstream")
                    cfg.bnd_condition.top_type = outstream;

                cfg.bnd_condition.top.x = node.attribute("u").as_double();
                cfg.bnd_condition.top.y = node.attribute("v").as_double();
                cfg.bnd_condition.top.z = node.attribute("w").as_double();
            }

            if (bc_node.child("Bottom") != NULL)
            {
                auto node = bc_node.child("Bottom");
                std::string type = node.attribute("type").value();

                if (type == "noslip")
                    cfg.bnd_condition.bottom_type = noslip;
                else if (type == "slip")
                    cfg.bnd_condition.bottom_type = slip;
                else if (type == "instream")
                    cfg.bnd_condition.bottom_type = instream;
                else if (type == "outstream")
                    cfg.bnd_condition.bottom_type = outstream;

                cfg.bnd_condition.bottom.x = node.attribute("u").as_double();
                cfg.bnd_condition.bottom.y = node.attribute("v").as_double();
                cfg.bnd_condition.bottom.z = node.attribute("w").as_double();
            }

            if (bc_node.child("Left") != NULL)
            {
                auto node = bc_node.child("Left");
                std::string type = node.attribute("type").value();

                if (type == "noslip")
                    cfg.bnd_condition.left_type = noslip;
                else if (type == "slip")
                    cfg.bnd_condition.left_type = slip;
                else if (type == "instream")
                    cfg.bnd_condition.left_type = instream;
                else if (type == "outstream")
                    cfg.bnd_condition.left_type = outstream;

                cfg.bnd_condition.left.x = node.attribute("u").as_double();
                cfg.bnd_condition.left.y = node.attribute("v").as_double();
                cfg.bnd_condition.left.z = node.attribute("w").as_double();
            }

            if (bc_node.child("Right") != NULL)
            {
                auto node = bc_node.child("Right");

                std::cout << node.attribute("type").value()  << std::endl;
                std::string type = node.attribute("type").value();

                if (type == "noslip")
                    cfg.bnd_condition.right_type = noslip;
                else if (type == "slip")
                    cfg.bnd_condition.right_type = slip;
                else if (type == "instream")
                    cfg.bnd_condition.right_type = instream;
                else if (type == "outstream")
                    cfg.bnd_condition.right_type = outstream;

                cfg.bnd_condition.right.x = node.attribute("u").as_double();
                cfg.bnd_condition.right.y = node.attribute("v").as_double();
                cfg.bnd_condition.right.z = node.attribute("w").as_double();
            }

            if (bc_node.child("Front") != NULL)
            {
                auto node = bc_node.child("Front");
                std::string type = node.attribute("type").value();

                if (type == "noslip")
                    cfg.bnd_condition.front_type = noslip;
                else if (type == "slip")
                    cfg.bnd_condition.front_type = slip;
                else if (type == "instream")
                    cfg.bnd_condition.front_type = instream;
                else if (type == "outstream")
                    cfg.bnd_condition.front_type = outstream;

                cfg.bnd_condition.front.x = node.attribute("u").as_double();
                cfg.bnd_condition.front.y = node.attribute("v").as_double();
                cfg.bnd_condition.front.z = node.attribute("w").as_double();
            }

            if (bc_node.child("Back") != NULL)
            {
                auto node = bc_node.child("Back");
                std::string type = node.attribute("type").value();

                if (type == "noslip")
                    cfg.bnd_condition.back_type = noslip;
                else if (type == "slip")
                    cfg.bnd_condition.back_type = slip;
                else if (type == "instream")
                    cfg.bnd_condition.back_type = instream;
                else if (type == "outstream")
                    cfg.bnd_condition.back_type = outstream;

                cfg.bnd_condition.back.x = node.attribute("u").as_double();
                cfg.bnd_condition.back.y = node.attribute("v").as_double();
                cfg.bnd_condition.back.z = node.attribute("w").as_double();
            }
        }
        else
        {
            std::cout << "No boundary conditions set!" << std::endl;
            std::exit(1);
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
            std::ifstream file(
                config_node.child("gridFile").first_attribute().as_string());

            if (!file)
            {
                std::cerr << "Could not open grid file at "
                    << config_node.child("gridFile").first_attribute().as_string()
                    << "!" << std::endl;
                std::exit(1);
            }

            cfg.num_fluid_cells = 0;

            std::size_t flag_res_x = cfg.num_x_blocks * cfg.cells_x_per_block + 2;
            std::size_t flag_res_y = cfg.num_y_blocks * cfg.cells_y_per_block + 2;
            std::size_t flag_res_z = cfg.num_y_blocks * cfg.cells_y_per_block + 2;

            cfg.flag_grid.resize(flag_res_x * flag_res_y * flag_res_z);

            std::size_t idx = (rank % (cfg.num_localities_x * cfg.num_localities_y)) % cfg.num_localities_x;
            std::size_t idy = (rank % (cfg.num_localities_x * cfg.num_localities_y)) / cfg.num_localities_x;
            std::size_t idz = rank / (cfg.num_localities_x * cfg.num_localities_y);

            std::size_t start_i = idx * cfg.num_x_blocks * cfg.cells_x_per_block;
            std::size_t end_i = start_i + cfg.num_x_blocks * cfg.cells_x_per_block;

            std::size_t start_j = idy * cfg.num_y_blocks * cfg.cells_y_per_block;
            std::size_t end_j = start_j + cfg.num_y_blocks * cfg.cells_y_per_block;

            std::size_t start_k = idz * cfg.num_z_blocks * cfg.cells_z_per_block;
            std::size_t end_k = start_k + cfg.num_z_blocks * cfg.cells_z_per_block;

            std::size_t offset_x = 0;
            std::size_t offset_y = 0;
            std::size_t offset_z = 0;

            std::size_t insert_i = 0;
            std::size_t insert_j = (end_j - start_j) - 1;
            std::size_t insert_k = (end_k - start_k) - 1;

            std::size_t i, j, k;
            i = 0;
            j = cfg.j_max + 1;
            k = cfg.k_max + 1;

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

                    std::bitset<9> flag(std::stoi(cell_val));

                    if (flag.test(is_fluid))
                        cfg.num_fluid_cells++;

                    if (i >= start_i && i < end_i && j >= start_j && j < end_j && k >= start_k && k < end_k)
                    {
                        cfg.flag_grid[(1 + insert_k) * flag_res_x * flag_res_y + (insert_j + 1)* flag_res_x + insert_i + 1] = flag;
                        ++insert_i;
                    }

                    if (!iss.good())
                        break;
                    ++i;
                }


                if (j >= start_j && j < end_j)
                    --insert_j;

                if (j == 0)
                {
                    insert_i = 0;
                    insert_j = (end_j - start_j) - 1;;
                    j = cfg.j_max + 1;

                    if (k >= start_k && k < end_k)
                        --insert_k;

                    --k;
                } else
                    --j;
            }
        }
        else
        {
            std::size_t flag_res_x = cfg.num_x_blocks * cfg.cells_x_per_block + 2;
            std::size_t flag_res_y = cfg.num_y_blocks * cfg.cells_y_per_block + 2;

            cfg.flag_grid.resize(flag_res_x * flag_res_y);
            //std::cerr << "Error: geometry file not set!" << std::endl;
           // std::exit(1);
        }

        if(config_node.child("initialUVFile") != NULL)
        {
            cfg.with_initial_uv_grid = true;

            std::ifstream file(
                config_node.child("initialUVFile")
                    .first_attribute().as_string());

            if (!file)
            {
                std::cerr << "Could not open initial velocity file at "
                    << config_node.child("initialUVFile").first_attribute().as_string()
                    << "!" << std::endl;
                std::exit(1);
            }

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
