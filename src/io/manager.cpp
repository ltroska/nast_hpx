#include <iostream>

#include "manager.hpp"
#include "pugixml/pugixml.hpp"

namespace io {

//!Method reads the config file given in XML format
/*!
 * @param path Path to XML config file
 */
config manager::read_config_from_file(const char *path)
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

	if(config_node.child("deltaT") != NULL)
	{
        config.delta_t = config_node.child("deltaT").first_attribute().as_double();
	}

	return config;
}
}//namespace io
