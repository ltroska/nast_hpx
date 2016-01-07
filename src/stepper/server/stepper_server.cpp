#include <hpx/include/iostreams.hpp>
#include <cmath>

#include "stepper_server.hpp"
#include "stepper/stencils.hpp"
#include "io/vtk_writer.hpp"

typedef stepper::server::stepper_server stepper_component;
typedef hpx::components::component<stepper_component> stepper_server_type;

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_COMPONENT(stepper_server_type, stepper_component);

HPX_REGISTER_ACTION(stepper::server::stepper_server::do_work_action, stepper_server_do_work_action);
HPX_REGISTER_ACTION(stepper::server::stepper_server::setup_action, stepper_server_setup_action);
HPX_REGISTER_ACTION(stepper::server::stepper_server::set_velocity_action, stepper_server_set_velocity_action);
HPX_REGISTER_ACTION(stepper::server::stepper_server::set_pressure_action, stepper_server_set_pressure_action);
HPX_REGISTER_ACTION(stepper::server::stepper_server::compute_fg_action, stepper_server_compute_fg_action);
HPX_REGISTER_ACTION(stepper::server::stepper_server::set_rhs_action, stepper_server_set_rhs_action);
HPX_REGISTER_ACTION(stepper::server::stepper_server::update_velocities_action, stepper_server_update_velocities_action);
HPX_REGISTER_ACTION(stepper::server::stepper_server::sor_cycle_action, stepper_server_sor_cycle_action);


namespace stepper { namespace server {

uint get_neighbor_id(uint id, direction dir, uint num_localities)
{
    uint res_x, res_y;
    if (num_localities == 2)
    {
        res_x = 2;
        res_y = 1;
    }
    else
    {
        res_x = static_cast<uint>(sqrt(num_localities));
        res_y = res_x;
    }

    switch (dir)
    {
        case left:
            return ((id-1)/res_x == id/res_x && id-1 < num_localities) ? id-1 : num_localities;

        case right:
            return ((id+1)/res_x == id/res_x && id+1 < num_localities) ? id+1 : num_localities;

        case top:
            return ((id+res_x) < num_localities && (id+res_x)/res_x == id/res_x +1) ? id+res_x : num_localities;

        case bottom:
            return ((id-res_x) < num_localities && (id-res_x)/res_x == id/res_x -1) ? id-res_x : num_localities;

        case top_left:
            return ((id+res_x-1) < num_localities && (id+res_x-1)/res_x == id/res_x +1) ? id+res_x-1 : num_localities;

        case top_right:
            return ((id+res_x+1) < num_localities && (id+res_x+1)/res_x == id/res_x+1) ? id+res_x+1 : num_localities;

        case bottom_left:
            return ((id-res_x-1) < num_localities && (id-res_x-1)/res_x == id/res_x-1) ? id-res_x-1 : num_localities;

        case bottom_right:
            return ((id-res_x+1) < num_localities && (id-res_x+1)/res_x == id/res_x-1) ? id-res_x+1 : num_localities;
        default:
            return num_localities;
    }
}

//get type of stepper
stepper_type get_type(uint id, uint num_localities)
{
    uint res = static_cast<uint>(sqrt(num_localities));

    uint id_mod_res = id % res;
    uint id_div_res = id / res;

    if (id_mod_res == 0)
    {
        if (id_div_res == 0)
            return stepper_type::bottom_left_boundary;
        if (id_div_res == res-1)
            return stepper_type::top_left_boundary;
        return stepper_type::left_boundary;
    }

    if (id_mod_res == res-1)
    {
        if (id_div_res == 0)
            return stepper_type::bottom_right_boundary;
        if (id_div_res == res-1)
            return stepper_type::top_right_boundary;
        return stepper_type::right_boundary;
    }

    if (id_div_res == 0)
        return stepper_type::bottom_boundary;

    if (id_div_res == res-1)
        return stepper_type::top_boundary;

    return stepper_type::interior;
}

stepper_server::stepper_server(uint num_localities)
    : num_localities_(num_localities),
      locality_id_(hpx::get_locality_id()),
      left_(hpx::find_from_basename(stepper_basename, get_neighbor_id(hpx::get_locality_id(), left, num_localities))),
      right_(hpx::find_from_basename(stepper_basename, get_neighbor_id(hpx::get_locality_id(), right, num_localities))),
      top_(hpx::find_from_basename(stepper_basename, get_neighbor_id(hpx::get_locality_id(), top, num_localities))),
      bottom_(hpx::find_from_basename(stepper_basename, get_neighbor_id(hpx::get_locality_id(), bottom, num_localities))),
      top_left_(hpx::find_from_basename(stepper_basename, get_neighbor_id(hpx::get_locality_id(), top_left, num_localities))),
      top_right_(hpx::find_from_basename(stepper_basename, get_neighbor_id(hpx::get_locality_id(), top_right, num_localities))),
      bottom_left_(hpx::find_from_basename(stepper_basename, get_neighbor_id(hpx::get_locality_id(), bottom_left, num_localities))),
      bottom_right_(hpx::find_from_basename(stepper_basename, get_neighbor_id(hpx::get_locality_id(), bottom_right, num_localities))),
      type_(get_type(hpx::get_locality_id(), num_localities))
{
 /*   hpx::cout << "on " << hpx::get_locality_id() << " left " << get_neighbor_id(hpx::get_locality_id(), left, num_localities) << " right " << get_neighbor_id(hpx::get_locality_id(), right, num_localities)
    << " top " << get_neighbor_id(hpx::get_locality_id(), top, num_localities) << " bottom " << get_neighbor_id(hpx::get_locality_id(), bottom, num_localities)
    << " top left "<<get_neighbor_id(hpx::get_locality_id(), top_left, num_localities) << " top right " << get_neighbor_id(hpx::get_locality_id(), top_right, num_localities)
    << " bottom left " << get_neighbor_id(hpx::get_locality_id(), bottom_left, num_localities) << " bottom right " << get_neighbor_id(hpx::get_locality_id(), bottom_right, num_localities)<< hpx::endl << hpx::flush;
*/}

uint stepper_server::setup(uint i_max, uint j_max, RealType x_length, RealType y_length, uint num_partitions_x, uint num_partitions_y)
{
    //padding for neighbor locality values
    num_local_partitions_x_ = num_partitions_x + 2;
    num_local_partitions_y_ = num_partitions_y + 2;

    x_length_ = x_length;
    y_length_ = y_length;

    dx_ = 1./i_max;
    dy_ = 1./j_max;

    i_max_ = i_max;
    j_max_ = j_max;

    //how many steppers for each dimension (i.e. if num_localities_ == 4, we need 2x2 steppers).
    if (num_localities_ == 2)
    {
        res_x_ = 2;
        res_y_ = 1;
    }
    else
    {
        res_x_ = static_cast<uint>(sqrt(num_localities_));
        res_y_ = res_x_;
    }

    num_cells_x_ = (i_max / res_x_) / num_partitions_x;
    num_cells_y_ = (j_max / res_y_) / num_partitions_x;

    U.resize(num_local_partitions_x_);
    for (uint i = 0; i < num_local_partitions_x_; ++i)
        U[i].resize(num_local_partitions_y_);

    uint num_local_cells_x, num_local_cells_y;

    for (uint j = 0; j < num_local_partitions_y_; ++j)
    {
        for (uint i = 0; i < num_local_partitions_x_; ++i)
        {
            num_local_cells_x = num_cells_x_;
            num_local_cells_y = num_cells_y_;
            if (i == 0 || i == num_local_partitions_x_ - 1)
                num_local_cells_x = 1;
            if (j == 0 || j == num_local_partitions_y_ - 1)
                num_local_cells_y = 1;
            U[i][j] = grid::partition(hpx::find_here(), num_local_cells_x, num_local_cells_y,
                                    ( i == 0 || j == 0 || i == num_local_partitions_x_ -1 || j == num_local_partitions_y_ -1 ) ? hpx::get_locality_id() +100 : hpx::get_locality_id() + j*0.1 + i*0.01);
        }
    }

   if (get_neighbor_id(hpx::get_locality_id(), left, num_localities_) < num_localities_)
    {
         //   hpx::cout << "sending left on " << hpx::get_locality_id() << hpx::endl << hpx::flush;

        for (uint j = 1; j < num_local_partitions_y_-1; j++)
            send_to_neighbor(j, U[1][j], left);
    }

    if (get_neighbor_id(hpx::get_locality_id(), right, num_localities_) < num_localities_)
    {
       // hpx::cout << "sending right on " << hpx::get_locality_id() << hpx::endl << hpx::flush;
        for (uint j = 1; j < num_local_partitions_y_-1; j++)
            send_to_neighbor(j, U[num_local_partitions_x_-2][j], right);
    }

    if (get_neighbor_id(hpx::get_locality_id(), bottom, num_localities_) < num_localities_)
    {
         //   hpx::cout << "sending bottom on " << hpx::get_locality_id() << hpx::endl << hpx::flush;

        for (uint i = 1; i < num_local_partitions_x_-1; i++)
            send_to_neighbor(i, U[i][1], bottom);
    }

    if (get_neighbor_id(hpx::get_locality_id(), top, num_localities_) < num_localities_)
    {
          //  hpx::cout << "sending top on " << hpx::get_locality_id() << hpx::endl << hpx::flush;

        for (uint i = 1; i < num_local_partitions_x_-1; i++)
            send_to_neighbor(i, U[i][num_local_partitions_y_-2], top);
    }

    if (get_neighbor_id(hpx::get_locality_id(), top_left, num_localities_) < num_localities_)
    {
          //  hpx::cout << "sending topleft on " << hpx::get_locality_id() << hpx::endl << hpx::flush;

        send_to_neighbor(0, U[1][num_local_partitions_y_-2], top_left);
    }

    if (get_neighbor_id(hpx::get_locality_id(), top_right, num_localities_) < num_localities_)
    {
           // hpx::cout << "sending topright on " << hpx::get_locality_id() << hpx::endl << hpx::flush;

        send_to_neighbor(0, U[num_local_partitions_x_-2][num_local_partitions_y_-2], top_right);
    }

    if (get_neighbor_id(hpx::get_locality_id(), bottom_left, num_localities_) < num_localities_)
    {
          //  hpx::cout << "sending bottomleft on " << hpx::get_locality_id() << hpx::endl << hpx::flush;

        send_to_neighbor(0, U[1][1], bottom_left);
    }

    if (get_neighbor_id(hpx::get_locality_id(), bottom_right, num_localities_) < num_localities_)
    {
          //  hpx::cout << "sending bottomright on " << hpx::get_locality_id() << hpx::endl << hpx::flush;

        send_to_neighbor(0, U[1][num_local_partitions_y_-2], bottom_right);
    }

    hpx::cout << "stepper on " << locality_id_ << " with " << num_cells_x_ << "x" << num_cells_y_ << " cells, " << num_local_partitions_x_
    << "x"<< num_local_partitions_y_ << " partitions, dx=" << dx_ << " dy=" << dy_ << hpx::endl << hpx::flush;
}

uint stepper_server::do_work(uint num_local_partitions_x, uint num_local_partitions_y, uint num_cells_x, uint num_cells_y, RealType delta_x, RealType delta_y)
{
    write_vtk_files();
    return 0;
}

uint stepper_server::set_velocity()
{
    /*
    *@TODO rewrite this to action
    */
    //left boundary
    if(hpx::get_locality_id() % res_x_ == 0)
    {
        for (uint k = 1; k < num_local_partitions_y_-1; k++)
        {
            //u_0,j=0 j=1,...,jmax
            grid::partition_data pdata = U[0][k].get_data(grid::center_partition).get();
            for (uint j = 0; j < num_cells_y_; j++)
            {
                grid::cell& cell = pdata.get_cell_ref(0,j);
                cell.u = 0;
            }

            //v_0,j=-v1,j j=1,...,jmax
            grid::partition_data pdata2 = U[1][k].get_data(grid::center_partition).get();
            for (uint j = 0; j < num_cells_y_; j++)
            {
                grid::cell& cell = pdata.get_cell_ref(0,j);
                cell.v = -pdata2.get_cell(0,j).v;
            }
        }
    }

    //right
    if(hpx::get_locality_id() % res_x_ == res_x_-1)
    {
        for (uint k = 1; k < num_local_partitions_y_-1; k++)
        {
            //u_imax,j=0 j=1,...,jmax
            grid::partition_data pdata = U[num_local_partitions_x_-2][k].get_data(grid::center_partition).get();
            for (uint j = 0; j < num_cells_y_; j++)
            {
                grid::cell& cell = pdata.get_cell_ref(num_cells_x_-1,j);
                cell.u = 0;
            }

            //v_imax+1,j=-v_imax,j j=1,...,jmax
            grid::partition_data pdata2 = U[num_local_partitions_x_-1][k].get_data(grid::center_partition).get();
            for (uint j = 0; j < num_cells_y_; j++)
            {
                grid::cell& cell = pdata2.get_cell_ref(0,j);
                cell.v = -pdata.get_cell(0,j).v;
            }
        }
    }

    //bottom
    if(hpx::get_locality_id() / res_x_ == 0)
    {
        for (uint k = 1; k < num_local_partitions_x_-1; k++)
        {
            //v_i,0=0 i=1,..,imax
            grid::partition_data pdata = U[k][0].get_data(grid::center_partition).get();
            for (uint i = 0; i < num_cells_x_; i++)
            {
                grid::cell& cell = pdata.get_cell_ref(i,0);
                cell.v = 0;
            }

            //u_i,0=-u_i,1 i=1,...,imax
            grid::partition_data pdata2=U[k][1].get_data(grid::center_partition).get();
            for (uint i = 0; i < num_cells_x_; i++)
            {
                grid::cell& cell = pdata.get_cell_ref(i,0);
                cell.u = -pdata2.get_cell(i,0).u;
            }
        }
    }

    //top
    if(hpx::get_locality_id() / res_x_ == res_y_-1)
    {
        for (uint k = 1; k < num_local_partitions_x_-1; k++)
        {
            //v_i,jmax=0 i=1,..,imax
            grid::partition_data pdata = U[k][num_local_partitions_y_-2].get_data(grid::center_partition).get();
            for (uint i = 0; i < num_cells_x_; i++)
            {
                grid::cell& cell = pdata.get_cell_ref(i,num_cells_y_-1);
                cell.v = 0;
            }

            //u_i,jmax+1=-u_i,jmax i=1,...,imax
            grid::partition_data pdata2 = U[k][num_local_partitions_y_-1].get_data(grid::center_partition).get();
            for (uint i = 0; i < num_cells_x_; i++)
            {
                grid::cell& cell = pdata2.get_cell_ref(i,0);
                //cell.u = -pdata.get_cell(i,num_cells_y_-1).u;

                //driven cavitz
                cell.u = 2.0-pdata.get_cell(i,num_cells_y_-1).u;
            }
        }
    }

    return 1;
}


uint stepper_server::set_pressure()
{
    /*
    *@rewrite this to action
    */

     //left boundary
    if(hpx::get_locality_id() % res_x_ == 0)
    {
        for (uint k = 1; k < num_local_partitions_y_-1; k++)
        {
            grid::partition_data pdata = U[0][k].get_data(grid::center_partition).get();
            grid::partition_data pdata2 = U[1][k].get_data(grid::center_partition).get();

            for (uint j = 0; j < num_cells_y_; j++)
            {
                grid::cell& cell = pdata.get_cell_ref(0,j);
                cell.p = pdata2.get_cell(0,j).p;
            }
        }
    }

     //right boundary
    if(hpx::get_locality_id() % res_x_ == res_x_-1)
    {
        for (uint k = 1; k < num_local_partitions_y_-1; k++)
        {
            grid::partition_data pdata = U[num_local_partitions_x_-1][k].get_data(grid::center_partition).get();
            grid::partition_data pdata2 = U[num_local_partitions_x_-2][k].get_data(grid::center_partition).get();

            for (uint j = 0; j < num_cells_y_; j++)
            {
                grid::cell& cell = pdata.get_cell_ref(0,j);
                cell.p = pdata2.get_cell(num_cells_x_-1,j).p;
            }
        }
    }

    //bottom boundary
    if(hpx::get_locality_id() / res_x_ == 0)
    {
        for (uint k = 1; k < num_local_partitions_x_-1; k++)
        {
            grid::partition_data pdata = U[k][0].get_data(grid::center_partition).get();
            grid::partition_data pdata2 = U[k][1].get_data(grid::center_partition).get();

            for (uint i = 0; i < num_cells_x_; i++)
            {
                grid::cell& cell = pdata.get_cell_ref(i,0);
                cell.p = pdata2.get_cell(i,0).p;
            }
        }
    }

    //top boundary
    if(hpx::get_locality_id() / res_x_ == res_y_-1)
    {
        for (uint k = 1; k < num_local_partitions_x_-1; k++)
        {
            grid::partition_data pdata = U[k][num_local_partitions_y_-1].get_data(grid::center_partition).get();
            grid::partition_data pdata2 = U[k][num_local_partitions_y_-2].get_data(grid::center_partition).get();

            for (uint i = 0; i < num_cells_x_; i++)
            {
                grid::cell& cell = pdata.get_cell_ref(i,0);
                cell.p = pdata2.get_cell(i,num_cells_y_-1).p;
            }
        }
    }


    return 1;
}

uint stepper_server::compute_fg()
{
    std::vector<hpx::future<uint> > fut;

    do_compute_fg_action act;

    for(uint i = 0; i < num_local_partitions_x_-1; i++)
        for(uint j = 0; j < num_local_partitions_y_-1; j++)
        {
            if(i==0 && j == 0)
                continue;

            fut.push_back(hpx::async(act, get_id(), i, j));
        }

    hpx::wait_all(fut);
    return 1;
}

uint stepper_server::do_compute_fg(uint i , uint j)
{
    grid::partition center = U[i][j];
    grid::partition_data pdata_center = center.get_data(grid::center_partition).get();

    if (hpx::get_locality_id() % res_x_ == 0 && i == 0)
    {
        for (uint l = 0; l < num_cells_y_; l++)
        {
            grid::cell& cell = pdata_center.get_cell_ref(0,l);
            cell.f = cell.u;
        }

        return 1;
    }



    if (hpx::get_locality_id() / res_x_ == 0 && j == 0)
    {
        for (uint k = 0; k < num_cells_x_; k++)
        {
            grid::cell& cell = pdata_center.get_cell_ref(k,0);
            cell.g = cell.v;
        }

        return 1;
    }

    if(i == 0 || j == 0 )
        return 1;

    grid::cell right, left, top, bottom, bottomright, topleft;

    for(uint k = 0; k < num_cells_x_; k++)
    {
        for(uint l = 0; l < num_cells_y_; l++)
        {
            grid::cell& middle = pdata_center.get_cell_ref(k, l);

            //equation 20
            //left
            if (hpx::get_locality_id() % res_x_ == 0 && i == 0 && k == 0)
            {
                middle.f = middle.u;
            }

            //right
            if (hpx::get_locality_id() % res_x_ ==  res_x_-1 && i == num_local_partitions_x_-2 && k == num_cells_x_-1)
            {
                middle.f = middle.u;
            }

            //bottom
            if (hpx::get_locality_id() / res_x_ == 0 && j == 0 && l == 0)
            {
                middle.g = middle.v;
            }

            //top
            if (hpx::get_locality_id() / res_x_ == res_y_-1 && j == num_local_partitions_y_-2 && l == num_cells_y_-1)
            {
                middle.g = middle.v;
            }


            if(k+1 < num_cells_x_)
                right = pdata_center.get_cell(k+1,l);
            else
                right = U[i+1][j].get_data(grid::right_partition).get().get_cell(0,l);

            if(k > 0)
                left = pdata_center.get_cell(k-1,l);
            else
                left = U[i-1][j].get_data(grid::left_partition).get().get_cell(0,l);

            if(l+1 < num_cells_y_)
                top = pdata_center.get_cell(k,l+1);
            else
                top = U[i][j+1].get_data(grid::top_partition).get().get_cell(k,0);

            if(l > 0)
                bottom = pdata_center.get_cell(k,l-1);
            else
                bottom = U[i][j-1].get_data(grid::bottom_partition).get().get_cell(k,0);

            if(k+1 < num_cells_x_)
            {
                if(l > 0)
                    bottomright = pdata_center.get_cell(k+1,l-1);
                else
                    bottomright = U[i][j-1].get_data(grid::center_partition).get().get_cell(k+1,0);
            }
            else
            {
                if(l > 0)
                    bottomright = U[i+1][j].get_data(grid::center_partition).get().get_cell(0,l-1);
                else
                    bottomright = U[i+1][j-1].get_data(grid::center_partition).get().get_cell(0,0);
            }

            if(k > 0)
            {
                if (l+1 < num_cells_y_)
                    topleft = pdata_center.get_cell(k-1, l+1);
                else
                    topleft = U[i][j+1].get_data(grid::center_partition).get().get_cell(k-1, 0);
            }
            else
            {
                if(l+1 < num_cells_y_)
                    topleft = U[i-1][j].get_data(grid::center_partition).get().get_cell(0, l+1);
                else
                    topleft = U[i-1][j+1].get_data(grid::center_partition).get().get_cell(0,0);
            }

            middle.f = middle.u + dt_*(
                                   1/Re_ * (second_derivative_fwd_bkwd_x(right.u, middle.u, left.u, dx_)
                                            + second_derivative_fwd_bkwd_y(top.u, middle.u, bottom.u, dy_))
                                            - first_derivative_of_square_x(right.u, middle.u, left.u, dx_, alpha_)
                                            - first_derivative_of_product_y(right.v, middle.v, bottom.v, bottomright.v,
                                                                        bottom.u, middle.u, top.u, dx_, alpha_)
                                     );

            middle.g = middle.v + dt_*(
                                   1/Re_ * (second_derivative_fwd_bkwd_x(right.v, middle.v, left.v, dx_)
                                            + second_derivative_fwd_bkwd_y(top.v, middle.v, bottom.v, dy_)
                                            - first_derivative_of_product_x(left.u, middle.u, top.u, topleft.u,
                                                                            left.v, middle.v, right.v, dx_, alpha_)
                                            - first_derivative_of_square_y(top.v, middle.v, bottom.v, dy_, alpha_)
                                            )
                                   );

        }
    }

    return 1;
}

uint stepper_server::set_rhs()
{
    std::vector<hpx::future<uint> > fut;

    do_set_rhs_action act;

    for(uint i = 1; i < num_local_partitions_x_-1; i++)
        for(uint j = 1; j < num_local_partitions_y_-1; j++)
        {
            fut.push_back(hpx::async(act, get_id(), i, j));
        }

    hpx::wait_all(fut);
    return 1;
}

uint stepper_server::do_set_rhs(uint i, uint j)
{
    grid::partition center = U[i][j];
    grid::partition_data pdata = center.get_data(grid::center_partition).get();

    grid::cell left, bottom;


    for(uint k = 0; k < num_cells_x_; k++)
    {
        for(uint l = 0; l < num_cells_y_; l++)
        {
            grid::cell& middle = pdata.get_cell_ref(k, l);

            if(k > 0)
                left = pdata.get_cell(k-1,l);
            else
                left = U[i-1][j].get_data(grid::left_partition).get().get_cell(0,l);

            if(l > 0)
                bottom = pdata.get_cell(k,l-1);
            else
                bottom = U[i][j-1].get_data(grid::bottom_partition).get().get_cell(k,0);


            middle.rhs = 1./dt_*( (middle.f - left.f)/dx_ + (middle.g - bottom.g)/dy_ );
        }
    }

    return 1;
}

uint stepper_server::update_velocities()
{
    std::vector<hpx::future<uint> > fut;

    do_update_velocities_action act;

    for(uint i = 1; i < num_local_partitions_x_-1; i++)
        for(uint j = 1; j < num_local_partitions_y_-1; j++)
        {

            fut.push_back(hpx::async(act, get_id(), i, j));
        }

    hpx::wait_all(fut);
    return 1;
}

uint stepper_server::do_update_velocities(uint i, uint j)
{
    grid::partition center = U[i][j];
    grid::partition_data pdata = center.get_data(grid::center_partition).get();

    grid::cell right, top;


    for(uint k = 0; k < num_cells_x_; k++)
    {
        for(uint l = 0; l < num_cells_y_; l++)
        {
            grid::cell& middle = pdata.get_cell_ref(k, l);

            if(k+1 < num_cells_x_)
                right = pdata.get_cell(k+1,l);
            else
                right = U[i+1][j].get_data(grid::right_partition).get().get_cell(0,l);

            if(l+1 < num_cells_y_)
                top = pdata.get_cell(k,l+1);
            else
                top = U[i][j+1].get_data(grid::top_partition).get().get_cell(k,0);

            middle.u = middle.f - dt_/dx_ * (right.p - middle.p);
            middle.v = middle.g - dt_/dy_ * (top.p  - middle.p);
        }
    }

    return 1;
}

uint stepper_server::sor_cycle()
{
    std::vector<hpx::future<uint> > fut;

    do_sor_cycle_action act;

    for(uint i = 1; i < num_local_partitions_x_-1; i++)
        for(uint j = 1; j < num_local_partitions_y_-1; j++)
        {

            fut.push_back(hpx::async(act, get_id(), i, j, 1));
        }

    hpx::wait_all(fut);
    return 1;
}

uint stepper_server::do_sor_cycle(uint i, uint j, uint odd)
{
    grid::partition center = U[i][j];
    grid::partition_data pdata = center.get_data(grid::center_partition).get();

    grid::cell right, left, top, bottom;

    uint startk;
    for(uint l = 0; l < num_cells_y_; l++)
    {
        if(odd)
        {
            startk = (hpx::get_locality_id()/res_x_*num_local_partitions_y_+l+1) % 2;
        }
        else
        {
            startk = (hpx::get_locality_id()/res_x_+l) % 2;
        }

        for(uint k = startk; k < num_cells_x_; k+=2)
        {
            grid::cell& middle = pdata.get_cell_ref(k,l);
            hpx::cout << k << " " << l << " " << middle.p << " " << hpx::endl <<  hpx::flush;
        }
    }

    return 1;
}

void stepper_server::write_vtk_files()
{
    std::vector<std::vector<grid::partition_data> > space_part;
    space_part.resize(num_local_partitions_x_);

    for (uint i = 0; i < num_local_partitions_x_; ++i)
    {
        space_part[i].resize(num_local_partitions_y_);
        for (uint j = 0; j < num_local_partitions_y_; ++j)
        {
            grid::partition_data base = U[i][j].get_data(grid::partition_type::center_partition).get();
            space_part[i][j] = grid::partition_data(base.size_x(), base.size_y());
            for (uint k = 0; k < base.size(); ++k) {
                space_part[i][j][k] = base[k];
            }
        }
    }

    boost::shared_ptr<hpx::lcos::local::promise<int> > p =  boost::make_shared<hpx::lcos::local::promise<int> >();
    io::do_async_write(space_part, num_local_partitions_x_, num_local_partitions_y_, num_cells_x_, num_cells_y_,p);
}

}//namespace server
}//namespace stepper
