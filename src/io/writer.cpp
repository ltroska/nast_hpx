#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <limits>

#include "writer.hpp"

namespace nast_hpx { namespace io {

void writer::write_vtk(grid_type const& p_data, grid_type const& u_data, grid_type const& v_data, type_grid const& cell_types,
    std::size_t res_x, std::size_t res_y, std::size_t i_max, std::size_t j_max,
    Real dx, Real dy, std::size_t step, std::size_t loc)
{
    std::size_t num_localities = res_x * res_y;

    std::size_t cells_x = p_data.size_x_ - 2;
    std::size_t cells_y = p_data.size_y_ - 2;
    std::size_t partitions_x = 1;
    std::size_t partitions_y = 1;

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
                << i_max + 1 << " 1 " << j_max + 1
                << " 0 0\" GhostLevel=\"0\">" << std::endl
            << "<PPointData>" << std::endl
          //  << "<DataArray type=\"Float32\" Name=\"vorticity\" />" << std::endl
          //  << "<DataArray type=\"Float32\" Name=\"strom\" />" << std::endl
          //  << "<DataArray type=\"Float32\" Name=\"heat\" />" << std::endl
            << "</PPointData>" << std::endl
            << "<PCellData>" << std::endl
            << "<DataArray type=\"Float32\" Name=\"pressure\" />" << std::endl
            << "<DataArray type=\"Int32\" Name=\"obstacle\" />" << std::endl
          //  << "<DataArray type=\"Float32\" Name=\"temperature\" />" << std::endl
            << "<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"2\" />"
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
            os  << "<Piece Extent=\""
                       << static_cast<int>(cells_x*partitions_x*(loc%res_x))
                << " " << cells_x*partitions_x*((loc%res_x)+1)
                << " " << static_cast<int>(cells_y*partitions_y*(loc/res_x))
                << " " << cells_y*partitions_y*((loc/res_x)+1)
                << " 0 0\" Source=\"field_" << step
                    << "_locality_" << loc << ".vtr\"></Piece>" << std::endl;
        }

        os << "</PRectilinearGrid>" << std::endl
            << "</VTKFile>" << std::endl;

        fb.close();
    }

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

    int start_x, end_x, start_y, end_y;

    start_x = cells_x*partitions_x*(loc%res_x) - 1;
    end_x = cells_x*partitions_x*((loc%res_x)+1) + 1;
    start_y = cells_y*partitions_y*(loc/res_x) - 1;
    end_y = cells_y*partitions_y*((loc/res_x)+1) + 1;


    std::string coordinate_x;
    std::string coordinate_y;
    std::string coordinate_z = "0\n";

    for (int x = - 1; x <= end_x - start_x; x++)
        coordinate_x += std::to_string(dx*(start_x + x)) + " ";

    coordinate_x += "\n";

    for (int y = - 1; y <= end_y - start_y; y++)
        coordinate_y += std::to_string(dy*(start_y + y)) + " ";

    coordinate_y += "\n";

    std::stringstream p_stream;
    p_stream << std::setprecision(std::numeric_limits<Real>::digits10);

    std::stringstream obstacle_stream;

    std::stringstream uv_stream;
    uv_stream << std::setprecision(std::numeric_limits<Real>::digits10);

    std::stringstream vorticity_stream;
    vorticity_stream
        << std::setprecision(std::numeric_limits<Real>::digits10);

    std::stringstream strom_stream;
    strom_stream << std::setprecision(std::numeric_limits<Real>::digits10);

    std::stringstream heat_stream;
    heat_stream << std::setprecision(std::numeric_limits<Real>::digits10);

    std::stringstream temp_stream;
    temp_stream << std::setprecision(std::numeric_limits<Real>::digits10);

    vorticity_stream << "0\n";
    vorticity_stream << "0\n";
    strom_stream << "0\n";
    strom_stream << "0\n";
    heat_stream << "0\n";
    heat_stream << "0\n";
    p_stream << "0\n";
    obstacle_stream << "0\n";
    temp_stream << "0\n";
    uv_stream << "0 0\n";

    for (uint i = 0; i < cells_x * partitions_x; i++)
    {
        p_stream << "0\n";
        obstacle_stream << "0\n";
        temp_stream << "0\n";
        uv_stream << "0 0\n";
        vorticity_stream << "0\n";
        strom_stream << "0\n";
        heat_stream << "0\n";
    }

    p_stream << "0\n";
    obstacle_stream << "0\n";
    temp_stream << "0\n";
    uv_stream << "0 0\n";

    vorticity_stream << "0\n";
    strom_stream << "0\n";
    heat_stream << "0\n";

    vorticity_stream << "0\n";
    vorticity_stream << "0\n";
    strom_stream << "0\n";
    strom_stream << "0\n";
    heat_stream << "0\n";
    heat_stream << "0\n";

    for (uint i = 0; i < cells_x * partitions_x; i++)
    {
        vorticity_stream << "0\n";
        strom_stream << "0\n";
        heat_stream << "0\n";
    }

    vorticity_stream << "0\n";
    strom_stream << "0\n";
    heat_stream << "0\n";

    for (uint l = 0; l < partitions_y; ++l)
    {
        for (uint j = 1; j < cells_y + 1; ++j)
        {
            p_stream << "0\n";
            obstacle_stream << "0\n";
            uv_stream << "0 0\n";

            for (uint k = 0; k < partitions_x; ++k)
            {
                for (uint i = 1; i < cells_x + 1; ++i)
                {
                    if (cell_types(i, j).test(is_fluid))
                        obstacle_stream << "0\n";
                    else
                        obstacle_stream << "1\n";

                    if (cell_types(i, j).count() == 1 || cell_types(i, j).none())
                    {
                        p_stream << "0\n";
                        uv_stream << "0 0\n";
                    }
                    else
                    {
                    p_stream
                        << p_data(i, j)  << "\n";


                    if (i == 0)
                        uv_stream << "0 ";
                    else
                        uv_stream << (u_data(i, j) + u_data(i - 1, j)) / 2. << " ";

                    if (j == 0)
                        uv_stream << "0\n";
                    else
                        uv_stream << (v_data(i, j) + v_data(i, j - 1)) / 2. << "\n";
                    }


                }
            }

            p_stream << "0\n";
            obstacle_stream << "0\n";
            uv_stream << "0 0\n";

        }

    }

    p_stream << "0\n";
    obstacle_stream << "0\n";
    uv_stream << "0 0\n";

    for (uint i = 0; i < cells_x * partitions_x; i++)
    {
        p_stream << "0\n";
    obstacle_stream << "0\n";

    uv_stream << "0 0\n";

    }
    uv_stream << "0 0\n";
    obstacle_stream << "0\n";

    p_stream << "0\n";


    std::string pdatastring = p_stream.str();
    std::string uvdatastring = uv_stream.str();
    std::string vorticitystring = vorticity_stream.str();
    std::string stromstring = strom_stream.str();
    std::string heatstring = heat_stream.str();
    std::string tempstring = temp_stream.str();

    os  << std::setprecision(std::numeric_limits<Real>::digits10)
        << "<?xml version=\"1.0\"?>" << std::endl
        << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
            << std::endl
        << "<RectilinearGrid WholeExtent=\""
            << start_x << " " << end_x << " " << start_y << " " << end_y
            << " 0 0\">" << std::endl
        << "<Piece Extent=\""
            << start_x << " " << end_x << " " << start_y << " " << end_y
            << " 0 0\">" << std::endl
        << "<PointData>" << std::endl
       // << "<DataArray type=\"Float32\" Name=\"vorticity\">" << std::endl
        //<< vorticitystring << std::endl
        //<< "</DataArray>" << std::endl
      /*  << "<DataArray type=\"Float32\" Name=\"strom\">" << std::endl
        << stromstring << std::endl
        << "</DataArray>" << std::endl
        << "<DataArray type=\"Float32\" Name=\"heat\">" << std::endl
        << heatstring << std::endl
        << "</DataArray>" << std::endl*/
        << "</PointData>" << std::endl
        << "<CellData>" << std::endl
        << "<DataArray type=\"Float32\" Name=\"pressure\">" << std::endl
        << pdatastring << std::endl
        << "</DataArray>" << std::endl
        << "<DataArray type=\"Int32\" Name=\"obstacle\">" << std::endl
        << obstacle_stream.str() << std::endl
        << "</DataArray>" << std::endl
        /*<< "<DataArray type=\"Float32\" Name=\"temperature\">" << std::endl
        << tempstring << std::endl
        << "</DataArray>" << std::endl*/
        << "<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"2\">"
            << std::endl
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
