#pragma once
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

void do_write_vtk(std::vector<std::vector<grid::partition_data<scalar_cell> > > const& p_grid,
                    std::vector<std::vector<grid::partition_data<vector_cell> > > const& uv_grid,
                    std::vector<std::vector<grid::partition_data<scalar_cell> > > const& stream_grid,
                    std::vector<std::vector<grid::partition_data<scalar_cell> > > const& vorticity_grid,
                    std::vector<std::vector<grid::partition_data<scalar_cell> > > const& heat_grid,
                    std::vector<std::vector<grid::partition_data<scalar_cell> > > const& temp_grid,
                    RealType dx, RealType dy, uint step, uint i_max, uint j_max,
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


    uint actual_partitions_x = partitions_x;
    uint actual_partitions_y = partitions_y;

    partitions_x -= 2;
    partitions_y -= 2;

    //master file
    if (hpx::get_locality_id() == 0)
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
            << "<VTKFile type=\"PRectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl
            << "<PRectilinearGrid WholeExtent=\"1 " << i_max + 1 << " 1 " << j_max + 1 << " 0 0\" GhostLevel=\"1\">" << std::endl
            << "<PPointData>" << std::endl
            << "<DataArray type=\"Float32\" Name=\"vorticity\" />" << std::endl
            << "<DataArray type=\"Float32\" Name=\"strom\" />" << std::endl
            << "<DataArray type=\"Float32\" Name=\"heat\" />" << std::endl
            << "</PPointData>" << std::endl
            << "<PCellData>" << std::endl
            << "<DataArray type=\"Float32\" Name=\"pressure\" />" << std::endl
            << "<DataArray type=\"Float32\" Name=\"temperature\" />" << std::endl
            << "<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"2\" />" << std::endl
            << "</PCellData>" << std::endl
            << "<PCoordinates>" << std::endl
            << "<PDataArray type=\"Float32\" Name=\"X_COORDINATES\" NumberOfComponents=\"1\"/>" << std::endl
            << "<PDataArray type=\"Float32\" Name=\"Y_COORDINATES\" NumberOfComponents=\"1\"/>" << std::endl
            << "<PDataArray type=\"Float32\" Name=\"Z_COORDINATES\" NumberOfComponents=\"1\"/>" << std::endl
            << "</PCoordinates>" << std::endl;

        for (uint loc = 0; loc < hpx::find_all_localities().size(); loc++)
        {
            os  << "<Piece Extent=\"" << static_cast<int>(cells_x*partitions_x*(loc%res_x))
                << " " << cells_x*partitions_x*((loc%res_x)+1)
                << " " << static_cast<int>(cells_y*partitions_y*(loc/res_x))
                << " " << cells_y*partitions_y*((loc/res_x)+1)
                << " 0 0\" Source=\"field_" << step << "_locality_" << loc << ".vtr\"></Piece>" << std::endl;
        }

        os << "</PRectilinearGrid>" << std::endl
            << "</VTKFile>" << std::endl;

        fb.close();
    }

    //slave file
    uint loc = hpx::get_locality_id();

    uint left = (loc%res_x == 0) ? 1 : 0;
    uint right = (loc%res_x == res_x - 1) ? 1 : 0;
    uint bottom = (loc/res_x == 0) ? 1 : 0;
    uint top = (loc/res_x == res_y-1) ? 1 : 0;

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

    uint end_k, end_l;

    std::stringstream p_stream;
    p_stream << std::setprecision(std::numeric_limits<RealType>::digits10);

    std::stringstream uv_stream;
    uv_stream << std::setprecision(std::numeric_limits<RealType>::digits10);

    std::stringstream vorticity_stream;
    vorticity_stream << std::setprecision(std::numeric_limits<RealType>::digits10);

    std::stringstream strom_stream;
    strom_stream << std::setprecision(std::numeric_limits<RealType>::digits10);

    std::stringstream heat_stream;
    heat_stream << std::setprecision(std::numeric_limits<RealType>::digits10);

    std::stringstream temp_stream;
    temp_stream << std::setprecision(std::numeric_limits<RealType>::digits10);

    uint k, l, i, j;

    k = 0;
    l = 0;
    i = 0;
    j = 0;

    uint cnt = 1;

    vorticity_stream << "0\n";
    vorticity_stream << "0\n";
    strom_stream << "0\n";
    strom_stream << "0\n";
    heat_stream << "0\n";
    heat_stream << "0\n";
    p_stream << "0\n";
    temp_stream << "0\n";
    uv_stream << "0 0\n";

    for (uint i = 0; i < cells_x * partitions_x; i++)
    {
        p_stream << "0\n";
        temp_stream << "0\n";
        uv_stream << "0 0\n";
        vorticity_stream << "0\n";
        strom_stream << "0\n";
        heat_stream << "0\n";
    }

    p_stream << "0\n";
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



    for (uint j = 1; j < actual_partitions_y - 1; j++)
    {
        for (int row = 0; row < cells_y; row++)
        {
            p_stream << "0\n";
            temp_stream << "0\n";
            uv_stream << "0 0\n";
            vorticity_stream << "0\n";
            vorticity_stream << "0\n";
            strom_stream << "0\n";
            strom_stream << "0\n";
            heat_stream << "0\n";
            heat_stream << "0\n";

            for (uint i = 1; i < actual_partitions_x - 1; i++)
            {
                for (int col = 0; col < cells_x; col++)
                {
                    p_stream << p_grid[i][j].get_cell(col, row).value << "\n";
                    temp_stream << temp_grid[i][j].get_cell(col, row).value << "\n";

                    if (!( (left && i == 1 && col == 0) || (bottom && j == 1 && row == 0) ))
                    {
                        RealType u;
                        RealType v;

                        if (col != 0)
                            u = (uv_grid[i][j].get_cell(col, row).first + uv_grid[i][j].get_cell(col - 1, row).first) / 2.;
                        else
                            u = (uv_grid[i - 1][j].get_cell(cells_x - 1, row).first + uv_grid[i][j].get_cell(0, row).first) / 2.;

                        if (row != 0)
                            v = (uv_grid[i][j].get_cell(col, row).second + uv_grid[i][j].get_cell(col, row - 1).second) / 2.;
                        else
                            v = (uv_grid[i][j - 1].get_cell(col, cells_y - 1).second + uv_grid[i][j].get_cell(col, 0).second) / 2.;

                        uv_stream << u << " " << v << "\n";

                    }
                    else
                    {
                        uv_stream << "0 0\n";
                    }


                    if (!( (left && i == 1 && col == 0) || (right && i == actual_partitions_x - 2 && (col == cells_x - 1 || col == cells_x - 2) )
                            || (bottom && j == 1 && row == 0) || (top && j == actual_partitions_y - 2 && (row == cells_y - 1 || row == cells_y - 2) )
                         )
                        )
                    {
                        vorticity_stream << vorticity_grid[i][j].get_cell(col, row).value << "\n";
                    }
                    else
                    {
                        vorticity_stream << "0\n";
                    }

                    if (!( (right && i == actual_partitions_x - 2 && col == cells_x - 1) || (top && j == actual_partitions_y - 2 && row == cells_y - 1) ) )
                    {
                        strom_stream << stream_grid[i][j].get_cell(col, row).value << "\n";
                        heat_stream << heat_grid[i][j].get_cell(col, row).value << "\n";
                    }
                    else
                    {
                        strom_stream << "0\n";
                        heat_stream << "0\n";
                    }
                }
            }

            p_stream << "0\n";
            temp_stream << "0\n";
            uv_stream << "0 0\n";
            vorticity_stream << "0\n";
            strom_stream << "0\n";
            heat_stream << "0\n";
        }

    }

    vorticity_stream << "0\n";
    vorticity_stream << "0\n";
    strom_stream << "0\n";
    strom_stream << "0\n";
    heat_stream << "0\n";
    heat_stream << "0\n";
    p_stream << "0\n";
    temp_stream << "0\n";
    uv_stream << "0 0\n";

    for (uint i = 0; i < cells_x * partitions_x; i++)
    {
        p_stream << "0\n";
        temp_stream << "0\n";
        uv_stream << "0 0\n";
        vorticity_stream << "0\n";
        strom_stream << "0\n";
        heat_stream << "0\n";
    }

    p_stream << "0\n";
    temp_stream << "0\n";
    uv_stream << "0 0\n";

    vorticity_stream << "0\n";
    strom_stream << "0\n";
    heat_stream << "0\n";


    std::string pdatastring = p_stream.str();
    std::string uvdatastring = uv_stream.str();
    std::string vorticitystring = vorticity_stream.str();
    std::string stromstring = strom_stream.str();
    std::string heatstring = heat_stream.str();
    std::string tempstring = temp_stream.str();

    os  << std::setprecision(std::numeric_limits<RealType>::digits10) << "<?xml version=\"1.0\"?>" << std::endl
        << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl
        << "<RectilinearGrid WholeExtent=\"" << start_x << " " << end_x << " " << start_y << " " << end_y
                << " 0 0\">" << std::endl
        << "<Piece Extent=\"" << start_x << " " << end_x << " " << start_y << " " << end_y
                << " 0 0\">" << std::endl
        << "<PointData>" << std::endl
        << "<DataArray type=\"Float32\" Name=\"vorticity\">" << std::endl
        << vorticitystring << std::endl
        << "</DataArray>" << std::endl
        << "<DataArray type=\"Float32\" Name=\"strom\">" << std::endl
        << stromstring << std::endl
        << "</DataArray>" << std::endl
        << "<DataArray type=\"Float32\" Name=\"heat\">" << std::endl
        << heatstring << std::endl
        << "</DataArray>" << std::endl
        << "</PointData>" << std::endl
        << "<CellData>" << std::endl
        << "<DataArray type=\"Float32\" Name=\"pressure\">" << std::endl
        << pdatastring << std::endl
        << "</DataArray>" << std::endl
        << "<DataArray type=\"Float32\" Name=\"temperature\">" << std::endl
        << tempstring << std::endl
        << "</DataArray>" << std::endl
        << "<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"2\">" << std::endl
        << uvdatastring << std::endl
        << "</DataArray>" << std::endl
        << "</CellData>" << std::endl
        << "<Coordinates>" << std::endl
        << "<DataArray type=\"Float32\" Name=\"X_COORDINATES\"  NumberOfComponents=\"1\" format=\"ascii\">" << std::endl
        << coordinate_x << std::endl
        << "</DataArray>" << std::endl
        << "<DataArray type=\"Float32\" Name=\"Y_COORDINATES\"  NumberOfComponents=\"1\" format=\"ascii\">" << std::endl
        << coordinate_y << std::endl
        << "</DataArray>" << std::endl
        << "<DataArray type=\"Float32\" Name=\"Z_COORDINATES\"  NumberOfComponents=\"1\" format=\"ascii\">" << std::endl
        << coordinate_z << std::endl
        << "</DataArray>" << std::endl
        << "</Coordinates>" << std::endl
        << "</Piece>" << std::endl
        << "</RectilinearGrid>" << std::endl
        << "</VTKFile>" << std::endl;

    fb.close();

    p->set_value(0);
}

// This function will be executed by an HPX thread
hpx::lcos::future<int> write_vtk_worker(std::vector<std::vector<grid::partition_data<scalar_cell> > > const& p_grid,
                    std::vector<std::vector<grid::partition_data<vector_cell> > > const& uv_grid,
                    std::vector<std::vector<grid::partition_data<scalar_cell> > > const& stream_grid,
                    std::vector<std::vector<grid::partition_data<scalar_cell> > > const& vorticity_grid,
                    std::vector<std::vector<grid::partition_data<scalar_cell> > > const& heat_grid,
                    std::vector<std::vector<grid::partition_data<scalar_cell> > > const& temp_grid,
                    RealType dx, RealType dy, uint step, uint i_max, uint j_max,
                        uint partitions_x, uint partitions_y, uint cells_x, uint cells_y)
{
    boost::shared_ptr<hpx::lcos::local::promise<int> > p =
        boost::make_shared<hpx::lcos::local::promise<int> >();

    // Get a reference to one of the IO specific HPX io_service objects ...
    hpx::threads::executors::io_pool_executor scheduler;

    // ... and schedule the handler to run on one of its OS-threads.
    scheduler.add(hpx::util::bind(&do_write_vtk, p_grid, uv_grid, stream_grid, vorticity_grid, heat_grid, temp_grid, dx, dy, step, i_max, j_max, partitions_x, partitions_y, cells_x, cells_y, p));

    return p->get_future();
}

int write_vtk(std::vector<std::vector<grid::partition_data<scalar_cell> > > const& p_grid,
                    std::vector<std::vector<grid::partition_data<vector_cell> > > const& uv_grid,
                    std::vector<std::vector<grid::partition_data<scalar_cell> > > const& stream_grid,
                    std::vector<std::vector<grid::partition_data<scalar_cell> > > const& vorticity_grid,
                    std::vector<std::vector<grid::partition_data<scalar_cell> > > const& heat_grid,
                    std::vector<std::vector<grid::partition_data<scalar_cell> > > const& temp_grid,
                    RealType dx, RealType dy, uint step, uint i_max, uint j_max,
                        uint partitions_x, uint partitions_y, uint cells_x, uint cells_y)
{
    hpx::lcos::future<int> f = write_vtk_worker(p_grid, uv_grid, stream_grid, vorticity_grid, heat_grid, temp_grid, dx, dy, step, i_max, j_max, partitions_x, partitions_y, cells_x, cells_y);
    return f.get();
}

}//namespace io

HPX_PLAIN_ACTION(io::write_vtk, write_vtk_action);


#endif
