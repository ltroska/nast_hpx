#include <hpx/lcos/gather.hpp>

#include <hpx/lcos/when_all.hpp>
#include <hpx/lcos/wait_all.hpp>
#include <hpx/lcos/when_each.hpp>

#include "partition_server.hpp"
#include "grid/stencils.hpp"
#include "io/writer.hpp"

typedef nast_hpx::grid::server::partition_server partition_component;
typedef hpx::components::component<partition_component> partition_server_type;

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_COMPONENT(partition_server_type, partition_component);

HPX_REGISTER_ACTION(nast_hpx::grid::server::partition_server::do_timestep_action,
    partition_server_do_timestep_action);
HPX_REGISTER_ACTION(nast_hpx::grid::server::partition_server::init_action,
    partition_server_init_action);
    
HPX_REGISTER_GATHER(Real, partition_server_residual_gather);
    
namespace nast_hpx { namespace grid { namespace server {
 
partition_server::partition_server(io::config&& cfg, uint idx, uint idy)
:   rank_(cfg.rank),
    num_partitions_(cfg.num_localities),
    num_partitions_x_(cfg.num_localities_x),
    num_partitions_y_(cfg.num_localities_y),
    cells_x_(cfg.cells_x_per_partition * cfg.num_local_partitions_x + 2),
    cells_y_(cfg.cells_y_per_partition * cfg.num_local_partitions_y + 2),
    idx_(idx),
    idy_(idy),
    is_left_(idx_ == 0),
    is_right_(idx_ == num_partitions_x_ - 1),
    is_bottom_(idy_ == 0),
    is_top_(idy_ == num_partitions_y_ - 1),
    num_x_blocks_(cfg.num_local_partitions_x),
    num_y_blocks_(cfg.num_local_partitions_y),
    cells_per_x_block_((cells_x_ - 2) / num_x_blocks_),
    cells_per_y_block_((cells_y_ - 2)/ num_y_blocks_),
    u_bnd_(cfg.u_bnd),
    v_bnd_(cfg.v_bnd),
    uv_bnd_type_(cfg.data_type),
    dx_(cfg.x_length / cfg.i_max),
    dy_(cfg.x_length / cfg.j_max),
    re_(cfg.re),
    alpha_(cfg.alpha),
    beta_(cfg.beta),
    eps_sq_(cfg.eps_sq),
    gx_(cfg.gx),
    gy_(cfg.gy),
    step_(0),
    over_dx_(1. / dx_),
    over_dy_(1. / dy_),
    dx_sq_(std::pow(dx_, 2)),
    dy_sq_(std::pow(dy_, 2)),
    over_dx_sq_(1. / dx_sq_),
    over_dy_sq_(1. / dy_sq_),
    part1_(1. - cfg.omega),
    part2_(cfg.omega * dx_sq_ * dy_sq_ / (2. * (dx_sq_ + dy_sq_))),
    num_fluid_cells_(cfg.num_fluid_cells),
    i_max_(cfg.i_max),
    j_max_(cfg.j_max),
    iter_max_(cfg.iter_max),
    vtk_(cfg.vtk),
    delta_vec_(cfg.delta_vec),
    next_out_(0)
{    
    std::cout << "dx " << dx_ << " dy " << dy_ << " re " << re_ << " gx " << gx_ << " gy " << gy_
        << " imax " << i_max_ << " jmax " << j_max_ << " alpha " << alpha_ << " beta " << beta_ << " part1 " << part1_
        << " ubnd " << u_bnd_ << " vbnd " << v_bnd_ << " type " << uv_bnd_type_ << " eps_sq " << eps_sq_ << std::endl;
    
    data_[U] = partition_data<Real>(cells_x_, cells_y_);
    data_[V] = partition_data<Real>(cells_x_, cells_y_);
    data_[F] = partition_data<Real>(cells_x_, cells_y_);
    data_[G] = partition_data<Real>(cells_x_, cells_y_);
    data_[P] = partition_data<Real>(cells_x_, cells_y_);
    rhs_data_ = partition_data<Real>(cells_x_, cells_y_);
    
    cell_types_.resize(cells_x_, cells_y_);
    
    for (std::size_t i = 0; i != cell_types_.size_; ++i)
        cell_types_[i] = cfg.flag_grid[i];
        
    cfg.flag_grid.clear();
}   

void partition_server::init()
{
    std::vector<hpx::future<hpx::id_type > > parts =
        hpx::find_all_from_basename(partition_basename, num_partitions_);
            
    ids_ = hpx::when_all(parts).then(hpx::util::unwrapped2(
                [](std::vector<hpx::id_type>&& ids) -> std::vector<hpx::id_type>
                { return ids;})
            ).get();
        
    if (!is_left_)
    {
        send_buffer_left_.dest_ = ids_[idy_ * num_partitions_x_ + idx_ - 1];
        
        recv_buffer_left_[P].valid_ = true;
        recv_buffer_left_[F].valid_ = true;
        recv_buffer_left_[G].valid_ = true;
        recv_buffer_left_[U].valid_ = true;
        recv_buffer_left_[V].valid_ = true;
        
        
        if (!is_top_)
        {
            send_buffer_top_left_.dest_ = ids_[(idy_ + 1) * num_partitions_x_ + idx_ - 1];
            
            recv_buffer_top_left_[U].valid_ = true;
            recv_buffer_top_left_[V].valid_ = true;            
        }
    }
    
    if (!is_right_)
    {
        send_buffer_right_.dest_ = ids_[idy_ * num_partitions_x_ + idx_ + 1];
        
        recv_buffer_right_[P].valid_ = true;
        recv_buffer_right_[U].valid_ = true;
        recv_buffer_right_[V].valid_ = true;
        
        if (!is_bottom_)
        {
            send_buffer_bottom_right_.dest_ = ids_[(idy_ - 1) * num_partitions_x_ + idx_ + 1];
            recv_buffer_bottom_right_[U].valid_ = true;
            recv_buffer_bottom_right_[V].valid_ = true;
        }
    }
    
    if (!is_bottom_)
    {
        send_buffer_bottom_.dest_ = ids_[(idy_ - 1) * num_partitions_x_ + idx_];
        
        recv_buffer_bottom_[P].valid_= true;
        recv_buffer_bottom_[F].valid_= true;
        recv_buffer_bottom_[G].valid_= true;
        recv_buffer_bottom_[U].valid_= true;
        recv_buffer_bottom_[V].valid_= true;
    }
    
    if (!is_top_)
    {
        send_buffer_top_.dest_ = ids_[(idy_ + 1) * num_partitions_x_ + idx_];
        
        recv_buffer_top_[P].valid_ = true;
        recv_buffer_top_[U].valid_ = true;
        recv_buffer_top_[V].valid_ = true;
    }
    
    for (uint i = 0; i < data_[P].size_x_ ; i++)
        for (uint j = 0; j < data_[P].size_y_; ++j)            
        {
          //  data_[U](i, j) = 100*(idy_*num_partitions_x_ + idx_)+i+0.01*j;
          /*  data_[V](i, j) = 100*(idy_*num_partitions_x_ + idx_)+i+0.01*j;
            data_[F](i, j) = 100*(idy_*num_partitions_x_ + idx_)+i+0.01*j;
            data_[G](i, j) = 100*(idy_*num_partitions_x_ + idx_)+i+0.01*j;
            data_[P](i, j) = 100*(idy_*num_partitions_x_ + idx_)+i+0.01*j;*/
        }
} 

template<>
void partition_server::send_boundary<LEFT>(std::size_t step, std::size_t var, future_vector& send_futures)
{
    if (!send_futures.empty())
        hpx::when_all(send_futures).then(
            hpx::launch::async,
            hpx::util::bind(
                boost::ref(send_buffer_left_),
                boost::ref(data_[var]),
                step,
                var
            )
        );
}

template<>
void partition_server::send_boundary<RIGHT>(std::size_t step, std::size_t var, future_vector& send_futures)
{
    if (!send_futures.empty())
        hpx::when_all(send_futures).then(
            hpx::launch::async,
            hpx::util::bind(
                boost::ref(send_buffer_right_),
                boost::ref(data_[var]),
                step,
                var
            )
        );
}

template<>
void partition_server::send_boundary<BOTTOM>(std::size_t step, std::size_t var, future_vector& send_futures)
{
    if (!send_futures.empty())
        hpx::when_all(send_futures).then(
            hpx::launch::async,
            hpx::util::bind(
                boost::ref(send_buffer_bottom_),
                boost::ref(data_[var]),
                step,
                var
            )
        );
}

template<>
void partition_server::send_boundary<TOP>(std::size_t step, std::size_t var, future_vector& send_futures)
{
    if (!send_futures.empty())
        hpx::when_all(send_futures).then(
            hpx::launch::async,
            hpx::util::bind(
                boost::ref(send_buffer_top_),
                boost::ref(data_[var]),
                step,
                var
            )
        );
}

template<>
void partition_server::send_boundary<BOTTOM_RIGHT>(std::size_t step, std::size_t var, future_vector& send_futures)
{
    if (!send_futures.empty())
        hpx::when_all(send_futures).then(
            hpx::launch::async,
            hpx::util::bind(
                boost::ref(send_buffer_bottom_right_),
                boost::ref(data_[var]),
                step,
                var
            )
        );
}

template<>
void partition_server::send_boundary<TOP_LEFT>(std::size_t step, std::size_t var, future_vector& send_futures)
{
    if (!send_futures.empty())
        hpx::when_all(send_futures).then(
            hpx::launch::async,
            hpx::util::bind(
                boost::ref(send_buffer_top_left_),
                boost::ref(data_[var]),
                step,
                var
            )
        );
}


template<std::size_t var>
void partition_server::send_right_and_top_boundaries(std::size_t step, future_grid& send_futures)
{
    send_boundary<RIGHT>(step, var, send_futures[RIGHT]);
    send_boundary<TOP>(step, var, send_futures[TOP]);
}

template<std::size_t var>
void partition_server::send_cross_boundaries(std::size_t step, future_grid& send_futures)
{
    send_boundary<LEFT>(step, var, send_futures[LEFT]);
    send_boundary<RIGHT>(step, var, send_futures[RIGHT]);
    send_boundary<BOTTOM>(step, var, send_futures[BOTTOM]);
    send_boundary<TOP>(step, var, send_futures[TOP]);
}

template<std::size_t var>
void partition_server::send_all_boundaries(std::size_t step, future_grid& send_futures)
{
    send_boundary<LEFT>(step, var, send_futures[LEFT]);
    send_boundary<RIGHT>(step, var, send_futures[RIGHT]);
    send_boundary<BOTTOM>(step, var, send_futures[BOTTOM]);
    send_boundary<TOP>(step, var, send_futures[TOP]);
    send_boundary<BOTTOM_RIGHT>(step, var, send_futures[BOTTOM_RIGHT]);
    send_boundary<TOP_LEFT>(step, var, send_futures[TOP_LEFT]);
}


template<>
void partition_server::receive_boundary<LEFT>(std::size_t step, std::size_t var, future_vector& recv_futures)
{
    if (recv_buffer_left_[var].valid_)   
        recv_futures[LEFT] =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_left_[var]),
                    boost::ref(data_[var]),
                    step
                )
        );
}

template<>
void partition_server::receive_boundary<RIGHT>(std::size_t step, std::size_t var, future_vector& recv_futures)
{
    if (recv_buffer_right_[var].valid_)   
        recv_futures[RIGHT] =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_right_[var]),
                    boost::ref(data_[var]),
                    step
                )
        );
}

template<>
void partition_server::receive_boundary<BOTTOM>(std::size_t step, std::size_t var, future_vector& recv_futures)
{
    if (recv_buffer_bottom_[var].valid_)   
        recv_futures[BOTTOM] =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_bottom_[var]),
                    boost::ref(data_[var]),
                    step
                )
        );
}

template<>
void partition_server::receive_boundary<TOP>(std::size_t step, std::size_t var, future_vector& recv_futures)
{
    if (recv_buffer_top_[var].valid_)   
        recv_futures[TOP] =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_top_[var]),
                    boost::ref(data_[var]),
                    step
                )
        );
}

template<>
void partition_server::receive_boundary<BOTTOM_RIGHT>(std::size_t step, std::size_t var, future_vector& recv_futures)
{
    if (recv_buffer_bottom_right_[var].valid_)   
        recv_futures[BOTTOM_RIGHT] =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_bottom_right_[var]),
                    boost::ref(data_[var]),
                    step
                )
        );
}

template<>
void partition_server::receive_boundary<TOP_LEFT>(std::size_t step, std::size_t var, future_vector& recv_futures)
{
    if (recv_buffer_top_left_[var].valid_)   
        recv_futures[TOP_LEFT] =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_top_left_[var]),
                    boost::ref(data_[var]),
                    step
                )
        );
}
    
template<std::size_t var>
void partition_server::receive_left_and_bottom_boundaries(std::size_t step, future_vector& recv_futures)
{
    receive_boundary<LEFT>(step, var, recv_futures);
    receive_boundary<BOTTOM>(step, var, recv_futures);
    
    recv_futures[RIGHT] = hpx::make_ready_future();
    recv_futures[TOP] = hpx::make_ready_future();
    recv_futures[BOTTOM_RIGHT] = hpx::make_ready_future();
    recv_futures[TOP_LEFT] = hpx::make_ready_future();
}  

template<std::size_t var>
void partition_server::receive_cross_boundaries(std::size_t step, future_vector& recv_futures)
{
    receive_boundary<LEFT>(step, var, recv_futures);
    receive_boundary<RIGHT>(step, var, recv_futures);
    receive_boundary<BOTTOM>(step, var, recv_futures);
    receive_boundary<TOP>(step, var, recv_futures);
    
    recv_futures[BOTTOM_RIGHT] = hpx::make_ready_future();
    recv_futures[TOP_LEFT] = hpx::make_ready_future();
}    

template<std::size_t var>
void partition_server::receive_all_boundaries(std::size_t step, future_vector& recv_futures)
{
    receive_boundary<LEFT>(step, var, recv_futures);
    receive_boundary<RIGHT>(step, var, recv_futures);
    receive_boundary<BOTTOM>(step, var, recv_futures);
    receive_boundary<TOP>(step, var, recv_futures);
    receive_boundary<BOTTOM_RIGHT>(step, var, recv_futures);
    receive_boundary<TOP_LEFT>(step, var, recv_futures);
}

void partition_server::wait_all_boundaries(future_vector& recv_futures)
{
  //  std::cout << "wait left" << std::endl;
    if (!is_left_)
    {
         //   std::cout << "wait left" << std::endl;

        recv_futures[LEFT].wait(); 
        
    }

    if (!is_right_)
    {
         //   std::cout << "wait RIGHT" << std::endl;

        recv_futures[RIGHT].wait(); 


    }
    if (!is_bottom_)
    {   // std::cout << "wait BOTTOM" << std::endl;

        recv_futures[BOTTOM].wait(); 
 
        if (!is_right_)
        {
          //      std::cout << "wait BOTTOM_RIGHT" << std::endl;

            recv_futures[BOTTOM_RIGHT].wait();
            
        }
    }
        
    if (!is_top_)
    {
          //  std::cout << "wait TOP" << std::endl;

        recv_futures[TOP].wait();

        if (!is_left_)
        {
             //   std::cout << "wait TOP_LEFT" << std::endl;

            recv_futures[TOP_LEFT].wait();
            
        }
    }
}

template<>
hpx::shared_future<void> partition_server::get_dependency<LEFT>(std::size_t idx_block,
    std::size_t idy_block, future_vector const& recv_futures,
    partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if (is_left_ && idx_block == 0)
        return hpx::make_ready_future();
    else if (idx_block == 0)
        return recv_futures[LEFT];
        
    return calc_futures(idx_block - 1, idy_block);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<RIGHT>(std::size_t idx_block,
    std::size_t idy_block, future_vector const& recv_futures,
    partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if (is_right_ && idx_block == calc_futures.size_x_ - 1)
        return hpx::make_ready_future();
    else if (idx_block == calc_futures.size_x_ - 1)
        return recv_futures[RIGHT];
        
    return calc_futures(idx_block + 1, idy_block);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<BOTTOM>(std::size_t idx_block,
    std::size_t idy_block, future_vector const& recv_futures,
    partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if (is_bottom_ && idy_block == 0)
        return hpx::make_ready_future();

    else if (idy_block == 0)
        return recv_futures[BOTTOM];
    
    return calc_futures(idx_block, idy_block - 1);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<TOP>(std::size_t idx_block,
    std::size_t idy_block, future_vector const& recv_futures,
    partition_data<hpx::shared_future<void> > const& calc_futures)
{   
    if (is_top_ && idy_block == calc_futures.size_y_ - 1)
        return hpx::make_ready_future();
    else if (idy_block == calc_futures.size_y_ - 1)
        return recv_futures[TOP];
                
    return calc_futures(idx_block, idy_block + 1);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<BOTTOM_RIGHT>(std::size_t idx_block,
    std::size_t idy_block, future_vector const& recv_futures,
    partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if ((is_bottom_ && idy_block == 0) || (is_right_ && idx_block == calc_futures.size_x_ - 1))
        return hpx::make_ready_future();
    
    if (idy_block == 0 && idx_block == calc_futures.size_x_ - 1)
        return recv_futures[BOTTOM_RIGHT];
    else if (idy_block == 0)
        return recv_futures[BOTTOM];
    else if (idx_block == calc_futures.size_x_ - 1)
        return recv_futures[RIGHT];
        
    return calc_futures(idx_block + 1, idy_block - 1);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<TOP_LEFT>(std::size_t idx_block,
    std::size_t idy_block, future_vector const& recv_futures,
    partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if ((is_left_ && idx_block == 0) || (is_top_ && idy_block == calc_futures.size_y_ - 1))
        return hpx::make_ready_future();
    
    if (idy_block == calc_futures.size_y_ - 1 && idx_block == 0)
        return recv_futures[TOP_LEFT];
    else if (idy_block == calc_futures.size_y_ - 1)
        return recv_futures[TOP];
    else if (idx_block == 0)
        return recv_futures[LEFT];
        
    return calc_futures(idx_block - 1, idy_block + 1);
}

std::pair<Real, Real> partition_server::do_timestep(Real dt)
{          
    
   // std::cout << "start step " << step_ << std::endl;
    hpx::util::high_resolution_timer t;
               
    partition_data<hpx::shared_future<void> > set_velocity_futures(num_x_blocks_, num_y_blocks_);
    
    future_grid send_futures_U(NUM_DIRECTIONS);
    future_grid send_futures_V(NUM_DIRECTIONS);
    
   // std::cout << "compute uv" << std::endl;            
  //  std::cout << "before U\n" << data_[U] << std::endl;
   // std::cout << "before V\n" << data_[V] << std::endl;      
  
    for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += cells_per_y_block_, ++ny_block)
    {
        range_type y_range(y, std::min(y + cells_per_y_block_, cells_y_ - 1));
    
        for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += cells_per_x_block_, ++nx_block)
        {
            range_type x_range(x, std::min(x + cells_per_x_block_, cells_x_ - 1));
                
           // std::cout << "block " << nx_block << " " << ny_block << " "
             //   << x_range.first << " " <<x_range.second <<" " <<y_range.first
               // << " " << y_range.second << std::endl;
            
            hpx::shared_future<void> calc_future =
                hpx::async(
                    hpx::util::bind(
                        &stencils<STENCIL_SET_VELOCITY>::call,
                        boost::ref(data_[U]), boost::ref(data_[V]),
                        boost::ref(data_[U]), boost::ref(data_[V]),
                        boost::ref(cell_types_), u_bnd_, v_bnd_,
                        x_range, y_range             
                    )
                );

            set_velocity_futures(nx_block, ny_block) = calc_future;

            if (!is_left_ && nx_block == 0)
            {
                send_futures_U[LEFT].push_back(calc_future);
                send_futures_V[LEFT].push_back(calc_future);
            }
                
            if (!is_right_ && nx_block == num_x_blocks_ - 1)
            {
                send_futures_U[RIGHT].push_back(calc_future);
                send_futures_V[RIGHT].push_back(calc_future);
            }

            if (!is_bottom_ && ny_block == 0)
            {
                send_futures_U[BOTTOM].push_back(calc_future);
                send_futures_V[BOTTOM].push_back(calc_future);
            }
        
            if (!is_top_ && ny_block == num_y_blocks_ - 1)
            {
                send_futures_U[TOP].push_back(calc_future);
                send_futures_V[TOP].push_back(calc_future);
            }
                
            if (!is_bottom_ && !is_right_ && nx_block == num_x_blocks_ - 1 && ny_block == 0)
            {
                send_futures_U[BOTTOM_RIGHT].push_back(calc_future);
                send_futures_V[BOTTOM_RIGHT].push_back(calc_future);
            }
                
            if (!is_top_ && !is_left_ && nx_block == 0 && ny_block == num_y_blocks_ - 1)
            {
                send_futures_U[TOP_LEFT].push_back(calc_future);
                send_futures_V[TOP_LEFT].push_back(calc_future);
            }              
        }
    }
    
    send_all_boundaries<U>(step_, send_futures_U);
    send_all_boundaries<V>(step_, send_futures_V);
    
    future_vector recv_futures_U(NUM_DIRECTIONS);
    future_vector recv_futures_V(NUM_DIRECTIONS);
                
    receive_all_boundaries<U>(step_, recv_futures_U);
    receive_all_boundaries<V>(step_, recv_futures_V);
    
    /*wait_all_boundaries(recv_futures_U);
    wait_all_boundaries(recv_futures_V);
    hpx::wait_all(set_velocity_futures.data_);
    
    std::cout << "after U\n" << data_[U] << std::endl;
    std::cout << "after V\n" << data_[V] << std::endl;*/

    if (vtk_ && next_out_ <= t_)
    {
        next_out_ += delta_vec_;
        
        hpx::when_all(set_velocity_futures.data_).then(
            hpx::launch::async,
            hpx::util::bind(
                &io::writer::write_vtk,
                boost::ref(data_[P]), boost::ref(data_[U]), boost::ref(data_[V]),
                num_partitions_x_, num_partitions_y_, i_max_, j_max_, dx_, dx_, outcount_++,
                rank_            
            )
        ).wait();        
    }
    
    partition_data<hpx::shared_future<void> > compute_fg_futures(num_x_blocks_, num_y_blocks_);
  //  std::cout << "compute fg" << std::endl;

    future_vector dependencies;
    future_grid send_futures_F(NUM_DIRECTIONS);
    future_grid send_futures_G(NUM_DIRECTIONS);

    for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += cells_per_y_block_, ++ny_block)
    {
        range_type y_range(y, std::min(y + cells_per_y_block_, cells_y_ - 1));
    
        for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += cells_per_x_block_, ++nx_block)
        {
            range_type x_range(x, std::min(x + cells_per_x_block_, cells_x_ - 1));
            
          //  std::cout << "block " << nx_block << " " << ny_block << " "
          //      << x_range.first << " " <<x_range.second <<" " <<y_range.first
          //      << " " << y_range.second << std::endl;
                
            dependencies.clear();
            dependencies.reserve(11);
            
            dependencies.push_back(set_velocity_futures(nx_block, ny_block));
            dependencies.push_back(get_dependency<LEFT>(nx_block, ny_block, recv_futures_U, set_velocity_futures));
            dependencies.push_back(get_dependency<LEFT>(nx_block, ny_block, recv_futures_V, set_velocity_futures));
            dependencies.push_back(get_dependency<RIGHT>(nx_block, ny_block, recv_futures_U, set_velocity_futures));
            dependencies.push_back(get_dependency<RIGHT>(nx_block, ny_block, recv_futures_V, set_velocity_futures));
            dependencies.push_back(get_dependency<BOTTOM>(nx_block, ny_block, recv_futures_U, set_velocity_futures));
            dependencies.push_back(get_dependency<BOTTOM>(nx_block, ny_block, recv_futures_V, set_velocity_futures));
            dependencies.push_back(get_dependency<TOP>(nx_block, ny_block, recv_futures_U, set_velocity_futures));
            dependencies.push_back(get_dependency<TOP>(nx_block, ny_block, recv_futures_V, set_velocity_futures));
            dependencies.push_back(get_dependency<TOP_LEFT>(nx_block, ny_block, recv_futures_U, set_velocity_futures));
            dependencies.push_back(get_dependency<BOTTOM_RIGHT>(nx_block, ny_block, recv_futures_V, set_velocity_futures));
            
            
            hpx::shared_future<void> calc_future =
                hpx::when_all(dependencies).then(
                    hpx::launch::async,
                    hpx::util::bind(
                        &stencils<STENCIL_COMPUTE_FG>::call,
                        boost::ref(data_[F]), boost::ref(data_[G]),
                        boost::ref(data_[U]), boost::ref(data_[V]),
                        boost::ref(cell_types_), re_, gx_, gy_, beta_, dx_,
                        dy_, dt, alpha_, x_range, y_range             
                    )
                );

            compute_fg_futures(nx_block, ny_block) = calc_future;

            if (!is_right_ && nx_block == num_x_blocks_ - 1)
            {
                send_futures_F[RIGHT].push_back(calc_future);
                send_futures_G[RIGHT].push_back(calc_future);
            }
        
            if (!is_top_ && ny_block == num_y_blocks_ - 1)
            {
                send_futures_F[TOP].push_back(calc_future);
                send_futures_G[TOP].push_back(calc_future);
            }
        }
    }
    
    send_right_and_top_boundaries<F>(step_, send_futures_F);
    send_right_and_top_boundaries<G>(step_, send_futures_G);
    
    future_vector recv_futures_F(NUM_DIRECTIONS);
    future_vector recv_futures_G(NUM_DIRECTIONS);
                
    receive_left_and_bottom_boundaries<F>(step_, recv_futures_F);
    receive_left_and_bottom_boundaries<G>(step_, recv_futures_G);
    
   /* wait_all_boundaries(recv_futures_F);
    wait_all_boundaries(recv_futures_G);
    
    hpx::wait_all(compute_fg_futures.data_);
    
    std::cout << "F\n" << data_[F] << std::endl;
    std::cout << "G\n" << data_[G] << std::endl;*/
    
        
    partition_data<hpx::shared_future<void> > compute_rhs_futures(num_x_blocks_, num_y_blocks_);
  //  std::cout << "compute rhs" << std::endl;

    for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += cells_per_y_block_, ++ny_block)
    {
        range_type y_range(y, std::min(y + cells_per_y_block_, cells_y_ - 1));
    
        for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += cells_per_x_block_, ++nx_block)
        {
            range_type x_range(x, std::min(x + cells_per_x_block_, cells_x_ - 1));
            
          // std::cout << "block " << nx_block << " " << ny_block << " "
           //     << x_range.first << " " <<x_range.second <<" " <<y_range.first
           //     << " " << y_range.second << std::endl;
                  
            dependencies.clear();
            dependencies.reserve(3);
            
            dependencies.push_back(compute_fg_futures(nx_block, ny_block));
            dependencies.push_back(get_dependency<LEFT>(nx_block, ny_block, recv_futures_F, compute_fg_futures));
            dependencies.push_back(get_dependency<BOTTOM>(nx_block, ny_block, recv_futures_G, compute_fg_futures));
                         
            compute_rhs_futures(nx_block, ny_block) =
                hpx::when_all(dependencies).then(
                    hpx::launch::async,
                    hpx::util::bind(
                        &stencils<STENCIL_COMPUTE_RHS>::call,
                        boost::ref(rhs_data_), boost::ref(data_[F]),
                        boost::ref(data_[G]), boost::ref(cell_types_),
                        dx_, dy_, dt, x_range, y_range             
                    )
                );
        }
    }    
       
   /*  hpx::wait_all(compute_rhs_futures.data_);
    
    std::cout << "RHS\n" << rhs_data_ << std::endl;*/
    
    //TODO make calc_futures members
    partition_data<hpx::shared_future<void> > set_p_futures(num_x_blocks_, num_y_blocks_);
    std::array<partition_data<hpx::shared_future<void> >, 2> sor_cycle_futures;
    
    std::array<future_vector, 2> recv_futures_P;
        
    hpx::util::high_resolution_timer t1;
            
    std::size_t current = 1;
    std::size_t last = 0;
            
    sor_cycle_futures[current].resize(num_x_blocks_, num_y_blocks_);
    sor_cycle_futures[last].resize(num_x_blocks_, num_y_blocks_);
    
    recv_futures_P[current].resize(NUM_DIRECTIONS);
    recv_futures_P[last].resize(NUM_DIRECTIONS);
            
    for (auto& a : sor_cycle_futures[last].data_)    
        a = hpx::make_ready_future();    
        
    for (auto& a : recv_futures_P[last])    
        a = hpx::make_ready_future();

        
    //std::cout << "start sor"<<std::endl;
        
    Real rres;
    std::size_t iter = 0;
    for (iter = 0; iter < iter_max_; ++iter)
    {
        std::vector<hpx::future<Real> > local_residuals;
        future_grid send_futures_P(NUM_DIRECTIONS);
        
       // hpx::wait_all(sor_cycle_futures[last].data_);
        //std::cout << "Before set P\n" << data_[P] << std::endl;

        
        receive_cross_boundaries<P>(iter, recv_futures_P[current]);
        
        for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += cells_per_y_block_, ++ny_block)
        {
            range_type y_range(y, std::min(y + cells_per_y_block_, cells_y_ - 1));
        
            for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += cells_per_x_block_, ++nx_block)
            {
                range_type x_range(x, std::min(x + cells_per_x_block_, cells_x_ - 1));
                
              //  std::cout << "block " << nx_block << " " << ny_block << " "
               //     << x_range.first << " " <<x_range.second <<" " <<y_range.first
              //      << " " << y_range.second << std::endl;
                                   
                dependencies.clear();

                dependencies.push_back(sor_cycle_futures[last](nx_block, ny_block));
               
                set_p_futures(nx_block, ny_block) = 
                    hpx::when_all(dependencies).then(
                        hpx::launch::async,
                        hpx::util::bind(
                            &stencils<STENCIL_SET_P>::call,
                            boost::ref(data_[P]), boost::ref(cell_types_),
                            nx_block, ny_block,
                            x_range, y_range             
                        )
                    );
            }
        }
        
      //  hpx::wait_all(set_p_futures.data_);
        
      //  std::cout << "After set P\n" << data_[P] << std::endl;
        
        for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += cells_per_y_block_, ++ny_block)
        {
            range_type y_range(y, std::min(y + cells_per_y_block_, cells_y_ - 1));
        
            for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += cells_per_x_block_, ++nx_block)
            {
                range_type x_range(x, std::min(x + cells_per_x_block_, cells_x_ - 1));
                
                dependencies.clear();
                dependencies.reserve(6);
                dependencies.push_back(set_p_futures(nx_block, ny_block));
                dependencies.push_back(compute_rhs_futures(nx_block, ny_block));
                dependencies.push_back(get_dependency<LEFT>(nx_block, ny_block, recv_futures_P[current], sor_cycle_futures[current]));
                dependencies.push_back(get_dependency<BOTTOM>(nx_block, ny_block, recv_futures_P[current], sor_cycle_futures[current]));
                dependencies.push_back(get_dependency<RIGHT>(nx_block, ny_block, recv_futures_P[last], sor_cycle_futures[last]));
                dependencies.push_back(get_dependency<TOP>(nx_block, ny_block, recv_futures_P[last], sor_cycle_futures[last]));

                hpx::shared_future<void> calc_future =
                    hpx::when_all(dependencies).then(
                        hpx::launch::async,
                        hpx::util::bind(
                            &stencils<STENCIL_SOR_CYCLE>::call,
                            boost::ref(data_[P]), boost::ref(rhs_data_),
                            boost::ref(cell_types_), part1_, part2_, dx_sq_, dy_sq_,
                            nx_block, ny_block, iter,
                            x_range, y_range             
                        )
                    );
                    
                sor_cycle_futures[current](nx_block, ny_block) = calc_future;
                    
                if (!is_left_ && nx_block == 0)
                    send_futures_P[LEFT].push_back(calc_future);                

                if (!is_right_ && nx_block == num_x_blocks_ - 1)
                    send_futures_P[RIGHT].push_back(calc_future);                

                if (!is_bottom_ && ny_block == 0)
                    send_futures_P[BOTTOM].push_back(calc_future);                

                if (!is_top_ && ny_block == num_y_blocks_ - 1)
                    send_futures_P[TOP].push_back(calc_future);                
            }
        }
        
       // hpx::wait_all(sor_cycle_futures[current].data_);
      //  std::cout << "after sor P\n" << data_[P] << std::endl;
        
        send_cross_boundaries<P>(iter, send_futures_P);
        
        for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += cells_per_y_block_, ++ny_block)
        {
            range_type y_range(y, std::min(y + cells_per_y_block_, cells_y_ - 1));
        
            for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += cells_per_x_block_, ++nx_block)
            {
                range_type x_range(x, std::min(x + cells_per_x_block_, cells_x_ - 1));
                
               // std::cout << "block " << nx_block << " " << ny_block << " "
               //   << x_range.first << " " <<x_range.second <<" " <<y_range.first
               //    << " " << y_range.second << std::endl;
                    
                dependencies.clear();
                dependencies.reserve(3);
                dependencies.push_back(sor_cycle_futures[current](nx_block, ny_block));
                dependencies.push_back(get_dependency<RIGHT>(nx_block, ny_block, recv_futures_P[current], sor_cycle_futures[current]));
                dependencies.push_back(get_dependency<TOP>(nx_block, ny_block, recv_futures_P[current], sor_cycle_futures[current]));
                
                
                local_residuals.push_back(
                    hpx::when_all(dependencies).then(
                        hpx::launch::async,
                        hpx::util::bind(
                            &stencils<STENCIL_SOR_RESIDUAL>::call,
                            boost::ref(data_[P]), boost::ref(rhs_data_),
                            boost::ref(cell_types_), over_dx_sq_, over_dy_sq_,
                            x_range, y_range             
                        )
                    )
                );
            }
        }
            
            
        std::swap(current, last);
            
        hpx::future<Real> local_residual =
            hpx::when_all(local_residuals).then(
                    hpx::util::unwrapped2(
                        [num_fluid_cells = num_fluid_cells_](std::vector<Real> residuals)
                        -> Real
                        {
                            Real sum = 0;

                            for (std::size_t i = 0; i < residuals.size(); ++i)
                                sum += residuals[i];

                            return sum / num_fluid_cells;
                        }
                    )
                );
                        
        if (hpx::get_locality_id() == 0)
        {
            hpx::future<std::vector<Real> > partial_residuals =
                hpx::lcos::gather_here(residual_basename,
                                        std::move(local_residual),
                                        num_partitions_, step_ * iter_max_ + iter);

            hpx::future<Real> residual =
                partial_residuals.then(
                    hpx::util::unwrapped(
                        [](std::vector<Real> local_residuals)
                            -> Real
                        {
                            Real residual = 0;

                            for (std::size_t i = 0; i < local_residuals.size(); ++i)
                                residual += local_residuals[i];

                            return residual;
                        }
                    )
                );
                
            if (iter == 249)
                rres = residual.get();

            // decide if SOR should keep running or not
            //hpx::lcos::broadcast_apply<set_keep_running_action>(
              //  localities, step * c.iter_max + iter,
                //(iter < c.iter_max && res > c.eps_sq));
        }
        // if not root locality, send residual to root locality
        else
            hpx::lcos::gather_there(residual_basename, std::move(local_residual),
                                        step_ * iter_max_ + iter).wait();
    }
    
    

  //  hpx::wait_all(sor_cycle_futures[last].data_);
   //  std::cout << "before2 U\n" << data_[U] << std::endl;
   // std::cout << "before2 V\n" << data_[V] << std::endl;
    Real elapsed = t1.elapsed();
        
    std::vector<hpx::future<std::pair<Real, Real> > > local_max_uvs;
    
    for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += cells_per_y_block_, ++ny_block)
    {
        range_type y_range(y, std::min(y + cells_per_y_block_, cells_y_ - 1));
    
        for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += cells_per_x_block_, ++nx_block)
        {
            range_type x_range(x, std::min(x + cells_per_x_block_, cells_x_ - 1));
            
           // std::cout << "block " << nx_block << " " << ny_block << " "
           //    << x_range.first << " " <<x_range.second <<" " <<y_range.first
           //    << " " << y_range.second << std::endl;
                
            dependencies.clear();
            dependencies.reserve(3);
            dependencies.push_back(get_dependency<RIGHT>(nx_block, ny_block, recv_futures_P[last], sor_cycle_futures[last]));
            dependencies.push_back(get_dependency<TOP>(nx_block, ny_block, recv_futures_P[last], sor_cycle_futures[last]));
            
            local_max_uvs.push_back(
                hpx::when_all(dependencies).then(
                    hpx::launch::async,
                    hpx::util::bind(
                        &stencils<STENCIL_UPDATE_VELOCITY>::call,
                        boost::ref(data_[U]), boost::ref(data_[V]), boost::ref(data_[F]),
                        boost::ref(data_[G]), boost::ref(data_[P]),
                        boost::ref(cell_types_), dt, over_dx_, over_dy_,
                        x_range, y_range             
                    )
                )
            );
        }
    }
    
    hpx::future<std::pair<Real, Real> > local_max_uv =
            hpx::when_all(local_max_uvs).then(
                    hpx::util::unwrapped2(
                        [](std::vector<std::pair<Real, Real> > max_uvs)
                        -> std::pair<Real, Real> 
                        {
                            std::pair<Real, Real>  max_uv(0, 0);

                            for (std::size_t i = 0; i < max_uvs.size(); ++i)
                            {
                                max_uv.first = max_uvs[i].first > max_uv.first ? max_uvs[i].first : max_uv.first;
                                max_uv.second = max_uvs[i].second > max_uv.second ? max_uvs[i].second : max_uv.second;
                            }
                                
                            return max_uv;
                        }
                    )
                );
            
    auto m = local_max_uv.get();
    
  //   std::cout << "after2 U\n" << data_[U] << std::endl;
  //  std::cout << "after2 V\n" << data_[V] << std::endl;
    
    if (hpx::get_locality_id() == 0)
        std::cout << "step " << step_ << " t " << t_ << " dt " << dt << " iter " << iter << " residual " << rres << " elapsed " << elapsed << " average " << elapsed / iter_max_ << std::endl;
        
    t_ += dt;
    ++step_;

    return m;    
}   
    
}
}
}