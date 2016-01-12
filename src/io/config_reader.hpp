#ifndef IO_CONFIG_READER_HPP
#define IO_CONFIG_READER_HPP

#include "internal/cfd_config.hpp"

namespace io {

class config_reader
{
public:
	//! Method to read the simulation parameters from an XML file and store it in the data object
	/*!
	 * @param path the path including the filename
	 * @return The struct with the simulation data filled with the data from the
	 *                 file.
	 */
	static cfd_config read_config_file(const char *path);

};

} //namespace IO
#endif
