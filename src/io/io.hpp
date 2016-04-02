#pragma once
#ifndef IO_IO_HPP
#define IO_IO_HPP

#include <bitset>
#include <sstream>
#include <fstream>

#include "pugixml/pugixml.hpp"

#include "config.hpp"

namespace io {

	//! Method to read the simulation parameters from an XML file and store it in the data object
	/*!
	 * @param path the path including the filename
	 * @return The struct with the simulation data filled with the data from the
	 *                 file.
	 */
	config read_config_from_file(const char *path)
    {
        pugi::xml_document doc;
        config config;

        pugi::xml_parse_result result = doc.load_file(path);

        if (!result) {
            std::cerr << "Error loading file: " << path << std::endl;
            exit(1);
        }

        pugi::xml_node config_node = doc.child("SimulationConfig");
        if (config_node == NULL) {
            std::cerr << "Error: A simulation configuration must be defined!" << std::endl;
            exit(1);
        }

        if(config_node.child("iMax") != NULL) {
            config.i_max = config_node.child("iMax").first_attribute().as_int();
        }
        else {
            std::cerr << "Error: iMax not set!" << std::endl;
            exit(1);
        }

        if(config_node.child("jMax") != NULL)
        {
            config.j_max = config_node.child("jMax").first_attribute().as_int();
        }
        else {
            std::cerr << "Error: jMax not set!" << std::endl;
            exit(1);
        }

        if(config_node.child("iRes") != NULL)
        {
            config.i_res = config_node.child("iRes").first_attribute().as_int();
        }
        else {
            std::cerr << "Error: iRes not set!" << std::endl;
            exit(1);
        }

        if(config_node.child("jRes") != NULL)
        {
            config.j_res = config_node.child("jRes").first_attribute().as_int();
        }
        else {
            std::cerr << "Error: jRes not set!" << std::endl;
            exit(1);
        }

        if(config_node.child("xLength") != NULL)
        {
            config.x_length = config_node.child("xLength").first_attribute().as_double();
        }
        else {
            std::cerr << "Error: xLength not set!" << std::endl;
            exit(1);
        }

        if(config_node.child("yLength") != NULL)
        {
            config.y_length = config_node.child("yLength").first_attribute().as_double();
        }
        else {
            std::cerr << "Error: yLength not set!" << std::endl;
            exit(1);
        }

        if(config_node.child("Re") != NULL)
        {
            config.re = config_node.child("Re").first_attribute().as_double();
        }

        if(config_node.child("omega") != NULL)
        {
            config.omega = config_node.child("omega").first_attribute().as_double();
        }

        if(config_node.child("tau") != NULL)
        {
            config.tau = config_node.child("tau").first_attribute().as_double();
        }

        if(config_node.child("eps") != NULL)
        {
            config.eps = config_node.child("eps").first_attribute().as_double();
            config.eps_sq = config.eps * config.eps;
        }

        if(config_node.child("alpha") != NULL)
        {
            config.alpha = config_node.child("alpha").first_attribute().as_double();
        }

        if(config_node.child("iterMax") != NULL)
        {
            config.iter_max = config_node.child("iterMax").first_attribute().as_int();
        }

        if(config_node.child("tEnd") != NULL)
        {
            config.t_end = config_node.child("tEnd").first_attribute().as_double();
        }

        if(config_node.child("dt") != NULL)
        {
            config.dt = config_node.child("dt").first_attribute().as_double();
        }
        else
            config.dt = 0.01;

        if(config_node.child("outputSkipSize") != NULL)
        {
            config.output_skip_size = config_node.child("outputSkipSize").first_attribute().as_int();
        }
        else
        {
            config.output_skip_size = 1;
        }

        if(config_node.child("subIterations") != NULL)
        {
            config.sub_iterations = config_node.child("subIterations").first_attribute().as_int();
        }
        else
        {
            config.sub_iterations = 1;
        }

        if(config_node.child("withForEach") != NULL)
        {
            config.wfe = config_node.child("withForEach").first_attribute().as_int();
        }
        else
        {
            config.wfe = 1;
        }

        if(config_node.child("vtk") != NULL)
        {
            config.vtk = (config_node.child("vtk").first_attribute().as_int() == 1);
        }
        else
        {
            config.vtk = 0;
        }

        if(config_node.child("GX") != NULL)
        {
            config.gx = config_node.child("GX").first_attribute().as_double();
        }
        else
        {
            config.gx = 0;
        }

        if(config_node.child("GY") != NULL)
        {
            config.gy = config_node.child("GY").first_attribute().as_double();
        }
        else
        {
            config.gy = 0;
        }

        if(config_node.child("UO") != NULL)
        {
            config.u_bnd.top = config_node.child("UO").first_attribute().as_double();
        }
        else
        {
            config.u_bnd.top = 0;
        }

        if(config_node.child("UB") != NULL)
        {
            config.u_bnd.bottom = config_node.child("UB").first_attribute().as_double();
        }
        else
        {
            config.u_bnd.bottom = 0;
        }

        if(config_node.child("UL") != NULL)
        {
            config.u_bnd.left = config_node.child("UL").first_attribute().as_double();
        }
        else
        {
            config.u_bnd.left = 0;
        }

        if(config_node.child("UR") != NULL)
        {
            config.u_bnd.right = config_node.child("UR").first_attribute().as_double();
        }
        else
        {
            config.u_bnd.right = 0;
        }

        if(config_node.child("VO") != NULL)
        {
            config.v_bnd.top = config_node.child("VO").first_attribute().as_double();
        }
        else
        {
            config.v_bnd.top = 0;
        }

        if(config_node.child("VB") != NULL)
        {
            config.v_bnd.bottom = config_node.child("VB").first_attribute().as_double();
        }
        else
        {
            config.v_bnd.bottom = 0;
        }

        if(config_node.child("VL") != NULL)
        {
            config.v_bnd.left = config_node.child("VL").first_attribute().as_double();
        }
        else
        {
            config.v_bnd.left = 0;
        }

        if(config_node.child("VR") != NULL)
        {
            config.v_bnd.right = config_node.child("VR").first_attribute().as_double();
        }
        else
        {
            config.v_bnd.right = 0;
        }

        if(config_node.child("TO") != NULL)
        {
            config.temp_bnd.top = config_node.child("TO").first_attribute().as_double();
        }
        else
        {
            config.temp_bnd.top = 0;
        }

        if(config_node.child("TB") != NULL)
        {
            config.temp_bnd.bottom = config_node.child("TB").first_attribute().as_double();
        }
        else
        {
            config.temp_bnd.bottom = 0;
        }

        if(config_node.child("TL") != NULL)
        {
            config.temp_bnd.left = config_node.child("TL").first_attribute().as_double();
        }
        else
        {
            config.temp_bnd.left = 0;
        }

        if(config_node.child("TR") != NULL)
        {
            config.temp_bnd.right = config_node.child("TR").first_attribute().as_double();
        }
        else
        {
            config.temp_bnd.right = 0;
        }

        if(config_node.child("WL") != NULL)
        {
            config.data_type.left = config_node.child("WL").first_attribute().as_double();
        }
        else
        {
            config.data_type.left = 1;
        }

        if(config_node.child("WR") != NULL)
        {
            config.data_type.right = config_node.child("WR").first_attribute().as_double();
        }
        else
        {
            config.data_type.right = 1;
        }

        if(config_node.child("WB") != NULL)
        {
            config.data_type.bottom = config_node.child("WB").first_attribute().as_double();
        }
        else
        {
            config.data_type.bottom = 1;
        }

        if(config_node.child("WT") != NULL)
        {
            config.data_type.top = config_node.child("WT").first_attribute().as_double();
        }
        else
        {
            config.data_type.top = 1;
        }

        if(config_node.child("deltaVec") != NULL)
        {
            config.delta_vec = config_node.child("deltaVec").first_attribute().as_double();
        }
        else
        {
            config.delta_vec = 0;
        }

        if(config_node.child("gridFile") != NULL)
        {
            config.with_flag_grid = true;

            std::ifstream file(config_node.child("gridFile").first_attribute().as_string());

            uint i_max = 0;
            uint j_max = 0;


            while (true)
            {
                std::string line;
                std::getline(file, line);

                if (!file.good())
                    break;
                i_max = 0;

                std::stringstream iss(line);

                while (true)
                {
                    std::string cell_val;
                    std::getline(iss, cell_val, ',');

                    config.flag_grid.push_back(std::bitset<5>(std::stoi(cell_val)));
                    i_max++;

                    if (!iss.good())
                        break;
                }

                j_max++;
            }

            config.i_max = i_max - 2;
            config.j_max = j_max - 2;

        }
        else
        {
            config.with_flag_grid = false;
        }

        return config;
    }

} //namespace IO

#endif
