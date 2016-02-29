#ifndef IO_MANAGER_HPP
#define IO_MANAGER_HPP

#include "config.hpp"

namespace io {

struct manager
{
  //7  typedef std::vector<std::vector<grid::partition_data<RealType> > > grid_type;

	//! Method to read the simulation parameters from an XML file and store it in the data object
	/*!
	 * @param path the path including the filename
	 * @return The struct with the simulation data filled with the data from the
	 *                 file.
	 */
	static config read_config_from_file(const char *path);

};

} //namespace IO

#endif
