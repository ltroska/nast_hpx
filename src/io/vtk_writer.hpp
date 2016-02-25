#ifndef IO_VTK_WRITER_HPP
#define IO_VTK_WRITER_HPP

#include <hpx/include/iostreams.hpp>
#include <hpx/include/runtime.hpp>
#include <hpx/include/util.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/thread_executors.hpp>

#include <fstream>
#include <sstream>

#include "stepper/server/stepper_server.hpp"

namespace io {

template<typename T>
void do_async_print(std::vector<std::vector<grid::partition_data<T> > > const& grid, std::string const message, uint partitions_x, uint partitions_y, uint cells_x, uint cells_y, boost::shared_ptr<hpx::lcos::local::promise<int> > p)
{
    uint res_x_, res_y_;
    if (hpx::get_num_localities_sync() == 2)
    {
        res_x_ = 2;
        res_y_ = 1;
    }
    else
    {
        res_x_ = static_cast<uint>(sqrt(hpx::get_num_localities_sync()));
        res_y_ = res_x_;
    }

    std::cout << message  << std::endl;
    for (uint j = partitions_y - 1; j <= partitions_y; --j)
    {
        for (uint row = cells_y - 1 ; row <= cells_y; --row)
        {
            for (uint i = 0; i < partitions_x; ++i)
            {
                for (uint col = 0; col <  cells_x; col++)
                {
                    std::cout << grid[i][j].get_cell(col, row) << " ";
                }
            }
            std::cout << std::endl;
        }
    }

    p->set_value(0);
}

void do_async_write(std::vector<std::vector<grid::partition_data<scalar_cell> > > const& p_grid,
                    std::vector<std::vector<grid::partition_data<vector_cell> > > const& uv_grid, RealType dx, RealType dy, uint step, uint i_max, uint j_max,
                        uint partitions_x, uint partitions_y, uint cells_x, uint cells_y, boost::shared_ptr<hpx::lcos::local::promise<int> > p)
{
    uint res_x, res_y;
    if (hpx::get_num_localities_sync() == 2)
    {
        res_x = 2;
        res_y = 1;
    }
    else
    {
        res_x = static_cast<uint>(sqrt(hpx::get_num_localities_sync()));
        res_y = res_x;
    }

    //master file
    if (hpx::get_locality_id() == 0)
    {
        std::string filename;
        filename.append ("./fields/");
        filename.append ("field_");
        filename.append (std::to_string(step));
        filename.append (".pvts");

        std::filebuf fb;
        fb.open (const_cast < char *>(filename.c_str ()), std::ios::out);
        std::ostream os (&fb);

        os  << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl
            << "<PStructuredGrid WholeExtent=\"0 " << i_max+1 << " 0 " << j_max+1 << " 0 0\" GhostLevel=\"1\">" << std::endl
            << "<PPointData Scalars=\"pressure\">" << std::endl
            << "<DataArray type=\"Float32\" Name=\"pressure\" />" << std::endl
            << "<DataArray type=\"Float32\" Name=\"vorticity\" />" << std::endl
            << "<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" />" << std::endl
            << "</PPointData>" << std::endl
            << "<PCellData></PCellData>" << std::endl
            << "<PPoints>" << std::endl
            << "<DataArray type=\"Float32\" Name=\"coordinates\"  NumberOfComponents=\"3\" format=\"ascii\" />" << std::endl
            << "</PPoints>" << std::endl;

        for (uint loc = 0; loc < hpx::find_all_localities().size(); loc++)
        {
            uint left = (loc%res_x == 0) ? 1 : 0;
            uint right = (loc%res_x == res_x - 1) ? 1 : 0;
            uint bottom = (loc/res_x == 0) ? 1 : 0;
            uint top = (loc/res_x == res_y-1) ? 1 : 0;

            os  << "<Piece Extent=\"" << cells_x*partitions_x*(loc%res_x)- (1-left)
                << " " << cells_x*partitions_x*((loc%res_x)+1)-right
                << " " << cells_y*partitions_y*(loc/res_x) - (1-bottom)
                << " " << cells_y*partitions_y*((loc/res_x)+1)-top
                << " 0 0\" Source=\"field_" << step << "_locality_" << loc << ".vts\"></Piece>" << std::endl;
        }

        os << "</PStructuredGrid>" << std::endl
            << "</VTKFile>" << std::endl;

        fb.close();
    }

    //slave file
    uint loc = hpx::get_locality_id();

    std::string filename;
    filename.append ("./fields/");
    filename.append ("field_");
    filename.append (std::to_string(step));
    filename.append ("_locality_");
    filename.append (std::to_string(loc));
    filename.append (".vts");

    std::filebuf fb;
    fb.open (const_cast < char *>(filename.c_str ()), std::ios::out);
    std::ostream os (&fb);

    uint left = (loc%res_x == 0) ? 1 : 0;
    uint right = (loc%res_x == res_x - 1) ? 1 : 0;
    uint bottom = (loc/res_x == 0) ? 1 : 0;
    uint top = (loc/res_x == res_y-1) ? 1 : 0;

    uint start_x, end_x, start_y, end_y;

    start_x = cells_x*partitions_x*(loc%res_x) - (1-left);
    end_x = cells_x*partitions_x*((loc%res_x)+1)-right;
    start_y = cells_y*partitions_y*(loc/res_x) - (1-bottom);
    end_y = cells_y*partitions_y*((loc/res_x)+1)-top;


    std::string coordinatestring;

    uint end_k, end_l;

    RealType p_data, u_data, v_data;

    std::stringstream p_stream;
    p_stream << std::setprecision(std::numeric_limits<RealType>::digits10);

    std::stringstream uv_stream;
    uv_stream << std::setprecision(std::numeric_limits<RealType>::digits10);

    std::stringstream vorticity_stream;
    vorticity_stream << std::setprecision(std::numeric_limits<RealType>::digits10);

    RealType part1, part2;

    for (uint i = 0; i < cells_x*partitions_x+2; i++)
    {
        p_stream << "0\n";
        uv_stream << "0 0 0\n";
        vorticity_stream << "0\n";
    }

    for (uint j = 0; j < partitions_y; j++)
    {
        for (int row = 0; row < cells_y; row++)
        {
            for (uint i = 0; i < partitions_x; i++)
            {
                for (int col = 0; col < cells_x; col++)
                {
                    if (row + 1 < cells_y)
                        p_stream << (p_grid[i][j].get_cell(col, row).value + p_grid[i][j].get_cell(col, row+1).value)/2 << "\n";
                    else if (j+1 < partitions_y)
                        p_stream <<  (p_grid[i][j].get_cell(col, row).value + p_grid[i][j+1].get_cell(col, 0).value)/2 << "\n";
                    else
                        p_stream << p_grid[i][j].get_cell(col, row).value << "\n";

                    if (row + 1 < cells_y)
                        part1 = uv_grid[i][j].get_cell(col, row+1).first - uv_grid[i][j].get_cell(col, row).first;
                    else if (j+1 < partitions_y)
                        part1 = uv_grid[i][j+1].get_cell(col, 0).first - uv_grid[i][j].get_cell(col, row).first;
                    else
                        part1 = 0;

                    if (col + 1 < cells_x)
                        part2 = uv_grid[i][j].get_cell(col+1, row).second - uv_grid[i][j].get_cell(col, row).second;
                    else if (i+1 < partitions_x)
                        part2 = uv_grid[i+1][j].get_cell(0, row).second - uv_grid[i][j].get_cell(col, row).second;
                    else
                        part2 = 0;

                    uv_stream << uv_grid[i][j].get_cell(col, row).first << " " << uv_grid[i][j].get_cell(col, row).second << " 0 \n";

                    vorticity_stream << part1/dy - part2/dx << "\n";

                    coordinatestring += std::to_string(dx*(start_x + i*cells_x + col)) + " ";

                    coordinatestring += std::to_string(dy*(start_y + j*cells_y + row)) + " ";

                    coordinatestring += "0\n";

                    if ( (col == cells_x - 1 && i == partitions_x - 1) || (row == cells_y -1 && j == partitions_y -1 ))
                    {
                        p_stream << "0\n";

                        uv_stream << "0 0 0\n";

                        vorticity_stream << "0\n";
                    }
                }

            }
        }
    }

    for (uint i = 0; i < cells_x*partitions_x+2; i++)
    {
        p_stream << "0\n";
        uv_stream << "0 0 0\n";
        vorticity_stream << "0\n";
    }

    std::string pdatastring = p_stream.str();
    std::string uvdatastring = uv_stream.str();
    std::string vorticitystring = vorticity_stream.str();

    os  << std::setprecision(std::numeric_limits<RealType>::digits10) << "<?xml version=\"1.0\"?>" << std::endl
        << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl
        << "<StructuredGrid WholeExtent=\"" << start_x << " " << end_x << " " << start_y << " " << end_y
                << " 0 0\">" << std::endl
        << "<Piece Extent=\"" << start_x << " " << end_x << " " << start_y << " " << end_y
                << " 0 0\">" << std::endl
        << "<PointData Scalars=\"pressure\">" << std::endl
        << "<DataArray type=\"Float32\" Name=\"pressure\">" << std::endl
        << pdatastring << std::endl
        << "</DataArray>" << std::endl
        << "<DataArray type=\"Float32\" Name=\"vorticity\">" << std::endl
        << vorticitystring << std::endl
        << "</DataArray>" << std::endl
        << "<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\">" << std::endl
        << uvdatastring << std::endl
        << "</DataArray>" << std::endl
        << "</PointData>" << std::endl
        << "<CellData></CellData>" << std::endl
        << "<Points>" << std::endl
        << "<DataArray type=\"Float32\" Name=\"coordinates\"  NumberOfComponents=\"3\" format=\"ascii\">" << std::endl
        << coordinatestring << std::endl
        << "</DataArray>" << std::endl
        << "</Points>" << std::endl
        << "</Piece>" << std::endl
        << "</StructuredGrid>" << std::endl
        << "</VTKFile>" << std::endl;

    fb.close();

    p->set_value(0);
}
}//namespace io

#endif
