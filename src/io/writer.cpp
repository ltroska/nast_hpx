#include "writer.hpp"

#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <limits>

namespace nast_hpx { namespace io {

void writer::write_vtk(grid_type const& p_data, grid_type const& u_data,
            grid_type const& v_data, grid_type const& w_data, type_grid const& cell_types,
            std::size_t res_x, std::size_t res_y, std::size_t res_z, std::size_t i_max,
            std::size_t j_max, std::size_t k_max, double dx, double dy, double dz, std::size_t step,
            std::size_t loc, std::size_t idx, std::size_t idy, std::size_t idz)
{
    std::size_t num_localities = res_x * res_y * res_z;

    std::size_t cells_x = p_data.size_x_ - 2;
    std::size_t cells_y = p_data.size_y_ - 2;
    std::size_t cells_z = p_data.size_z_ - 2;
    std::size_t partitions_x = 1;
    std::size_t partitions_y = 1;
    std::size_t partitions_z = 1;

    if (loc == 0)
    {
        std::string filename;
        filename.append ("./fields/");
        filename.append ("field_");
        filename.append (std::to_string(step));
        filename.append (".pvtr");

        std::filebuf fb;
        fb.open (const_cast < char *>(filename.c_str ()), std::ios::out);
        std::ostream os (&fb);

        os  << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"PRectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
                << std::endl
            << "<PRectilinearGrid WholeExtent=\"1 "
                << i_max + 1 << " 1 " << j_max + 1 << " 1 " << k_max + 1
                << "\" GhostLevel=\"0\">" << std::endl
            << "<PPointData>" << std::endl
          //  << "<DataArray type=\"Float32\" Name=\"vorticity\" />" << std::endl
          //  << "<DataArray type=\"Float32\" Name=\"strom\" />" << std::endl
          //  << "<DataArray type=\"Float32\" Name=\"heat\" />" << std::endl
            << "</PPointData>" << std::endl
            << "<PCellData>" << std::endl
            << "<DataArray type=\"Float32\" Name=\"pressure\" />" << std::endl
            << "<DataArray type=\"Int32\" Name=\"is_fluid\" />" << std::endl
          //  << "<DataArray type=\"Float32\" Name=\"temperature\" />" << std::endl
            << "<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" />"
            << std::endl
            << "</PCellData>" << std::endl
            << "<PCoordinates>" << std::endl
            << "<PDataArray type=\"Float32\" Name=\"X_COORDINATES\" NumberOfComponents=\"1\"/>"
                << std::endl
            << "<PDataArray type=\"Float32\" Name=\"Y_COORDINATES\" NumberOfComponents=\"1\"/>"
                << std::endl
            << "<PDataArray type=\"Float32\" Name=\"Z_COORDINATES\" NumberOfComponents=\"1\"/>"
                << std::endl
            << "</PCoordinates>" << std::endl;

        for (uint loc = 0; loc < num_localities; loc++)
        {
            std::size_t tmp_idx = (loc % (res_x * res_y)) % res_x;
            std::size_t tmp_idy = (loc % (res_x * res_y)) / res_x;
            std::size_t tmp_idz = loc / (res_x * res_y);

            os  << "<Piece Extent=\""
                       << static_cast<int>(cells_x * partitions_x * tmp_idx)
                << " " << cells_x * partitions_x * (tmp_idx + 1)
                << " " << static_cast<int>(cells_y * partitions_y * tmp_idy)
                << " " << cells_y * partitions_y * (tmp_idy + 1)
                << " " << static_cast<int>(cells_z * partitions_y * tmp_idz)
                << " " << cells_y * partitions_y * (tmp_idz + 1)
                << "\" Source=\"field_" << step
                    << "_locality_" << loc << ".vtr\"></Piece>" << std::endl;
        }

        os << "</PRectilinearGrid>" << std::endl
            << "</VTKFile>" << std::endl;

        fb.close();
    }

    std::vector<double> strom(cells_x + 2, 0);

    std::string filename;
    filename.append ("./fields/");
    filename.append ("field_");
    filename.append (std::to_string(step));
    filename.append ("_locality_");
    filename.append (std::to_string(loc));
    filename.append (".vtr");

    std::filebuf fb;
    fb.open (const_cast < char *>(filename.c_str ()), std::ios::out);
    std::ostream os (&fb);

    int start_x, end_x, start_y, end_y, start_z, end_z;

    start_x = cells_x * partitions_x * idx - 1;
    end_x = cells_x * partitions_x * (idx + 1) + 1;
    start_y = cells_y * partitions_y * idy - 1;
    end_y = cells_y * partitions_y * (idy + 1) + 1;
    start_z = cells_z * partitions_z * idz - 1;
    end_z = cells_z * partitions_z * (idz + 1) + 1;

    std::string coordinate_x;
    std::string coordinate_y;
    std::string coordinate_z;

    for (int x = - 1; x <= end_x - start_x; ++x)
        coordinate_x += std::to_string(dx * (start_x + x)) + " ";

    coordinate_x += "\n";

    for (int y = - 1; y <= end_y - start_y; ++y)
        coordinate_y += std::to_string(dy * (start_y + y)) + " ";

    coordinate_y += "\n";

    for (int z = - 1; z <= end_z - start_z; ++z)
        coordinate_z += std::to_string(dz * (start_z + z)) + " ";

    coordinate_z += "\n";

    std::stringstream p_stream;
    p_stream << std::setprecision(std::numeric_limits<double>::digits10);

    std::stringstream fluid_stream;

    std::stringstream uv_stream;
    uv_stream << std::setprecision(std::numeric_limits<double>::digits10);

    std::stringstream vorticity_stream;
    vorticity_stream
        << std::setprecision(std::numeric_limits<double>::digits10);

    std::stringstream strom_stream;
    strom_stream << std::setprecision(std::numeric_limits<double>::digits10);

    std::stringstream heat_stream;
    heat_stream << std::setprecision(std::numeric_limits<double>::digits10);

    std::stringstream temp_stream;
    temp_stream << std::setprecision(std::numeric_limits<double>::digits10);

    for (uint i = 0; i < (cells_x * partitions_x + 2) * (cells_y * partitions_y + 2); i++)
    {
        uv_stream << "0 0 0\n";
        p_stream << "0\n";
        fluid_stream << "0\n";
    }

    for (uint l = 0; l < partitions_z; ++l)
    {
        for (uint k = 1; k <= cells_z; ++k)
        {

            for (uint i = 0; i < cells_x * partitions_x + 2; i++)
            {
                uv_stream << "0 0 0\n";
                p_stream << "0\n";
                fluid_stream << "0\n";
            }


            for (uint m = 0; m < partitions_y; ++m)
            {
                for (uint j = 1; j <= cells_y; ++j)
                {
                    uv_stream << "0 0 0\n";
                    p_stream << "0\n";
                    fluid_stream << "0\n";


                    for (uint n = 0; n < partitions_x; ++n)
                    {
                        for (uint i = 1; i <= cells_x; ++i)
                        {
                           // uv_stream << u_data(i, j, k) << " " << v_data(i, j, k) << " " << w_data(i, j, k) << "\n";
                          //  p_stream << cell_types(i, j, k).to_ulong() << "\n";
                            if (cell_types(i, j, k).test(is_fluid))
                                fluid_stream << "1\n";
                            else
                                fluid_stream << "0\n";


                            if (cell_types(i, j, k).test(is_fluid))
                            {
                                p_stream << p_data(i, j, k)  << "\n";

                                uv_stream << (u_data(i, j, k) + u_data(i - 1, j, k)) / 2. << " ";

                                uv_stream << (v_data(i, j, k) + v_data(i, j - 1, k)) / 2. << " ";

                                uv_stream << (w_data(i, j, k) + w_data(i, j, k - 1)) / 2. << "\n";
                            }
                            else
                            {
                                p_stream << "0\n";
                                uv_stream << "0 0 0\n";
                            }


                    /*        if (cell_types(i, j, k).test(is_fluid))
                            {
                                double tmp = strom[i] + u_data(i, j, k) * dy;
                                strom_stream << tmp << "\n";

                                strom[i] = tmp;
                            }
                            else
                            {
                                strom_stream << "0\n";
                            }*/


                        }
                    }

                    uv_stream << "0 0 0\n";
                    p_stream << "0\n";
                    fluid_stream << "0\n";

                }
            }


            for (uint i = 0; i < cells_x * partitions_x + 2; i++)
            {
                uv_stream << "0 0 0\n";
                p_stream << "0\n";
                fluid_stream << "0\n";

            }
        }
    }

    for (uint i = 0; i < (cells_x * partitions_x + 2) * (cells_y * partitions_y + 2); i++)
    {
        uv_stream << "0 0 0\n";
        p_stream << "0\n";
        fluid_stream << "0\n";
    }


    std::string pdatastring = p_stream.str();
    std::string uvdatastring = uv_stream.str();
    std::string vorticitystring = vorticity_stream.str();
    std::string stromstring = strom_stream.str();
    std::string heatstring = heat_stream.str();
    std::string tempstring = temp_stream.str();

    os  << std::setprecision(std::numeric_limits<double>::digits10)
        << "<?xml version=\"1.0\"?>" << std::endl
        << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
            << std::endl
        << "<RectilinearGrid WholeExtent=\""
            << start_x << " " << end_x << " " << start_y << " " << end_y << " " << start_z << " " << end_z << "\">" << std::endl
        << "<Piece Extent=\""
            << start_x << " " << end_x << " " << start_y << " " << end_y << " " << start_z << " " << end_z << "\">" << std::endl
        << "<PointData>" << std::endl
       // << "<DataArray type=\"Float32\" Name=\"vorticity\">" << std::endl
        //<< vorticitystring << std::endl
        //<< "</DataArray>" << std::endl
      //  << "<DataArray type=\"Float32\" Name=\"strom\">" << std::endl
      //  << stromstring << std::endl
      //  << "</DataArray>" << std::endl/*
      //  << "<DataArray type=\"Float32\" Name=\"heat\">" << std::endl
      //  << heatstring << std::endl
       // << "</DataArray>" << std::endl*/
        << "</PointData>" << std::endl
        << "<CellData>" << std::endl
        << "<DataArray type=\"Float32\" Name=\"pressure\">" << std::endl
        << pdatastring << std::endl
        << "</DataArray>" << std::endl
        << "<DataArray type=\"Int32\" Name=\"is_fluid\">" << std::endl
        << fluid_stream.str() << std::endl
        << "</DataArray>" << std::endl
        /*<< "<DataArray type=\"Float32\" Name=\"temperature\">" << std::endl
        << tempstring << std::endl
        << "</DataArray>" << std::endl*/
        << "<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\">" << std::endl
        << uvdatastring << std::endl
        << "</DataArray>" << std::endl
        << "</CellData>" << std::endl
        << "<Coordinates>" << std::endl
        << "<DataArray type=\"Float32\" Name=\"X_COORDINATES\"  NumberOfComponents=\"1\" format=\"ascii\">"
            << std::endl
        << coordinate_x << std::endl
        << "</DataArray>" << std::endl
        << "<DataArray type=\"Float32\" Name=\"Y_COORDINATES\"  NumberOfComponents=\"1\" format=\"ascii\">"
            << std::endl
        << coordinate_y << std::endl
        << "</DataArray>" << std::endl
        << "<DataArray type=\"Float32\" Name=\"Z_COORDINATES\"  NumberOfComponents=\"1\" format=\"ascii\">"
            << std::endl
        << coordinate_z << std::endl
        << "</DataArray>" << std::endl
        << "</Coordinates>" << std::endl
        << "</Piece>" << std::endl
        << "</RectilinearGrid>" << std::endl
        << "</VTKFile>" << std::endl;

    fb.close();
}

}
}
