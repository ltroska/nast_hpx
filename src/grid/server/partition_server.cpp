

#include "partition_server.hpp"
#include "grid/stencils.hpp"

#include <hpx/lcos/when_all.hpp>
#include <hpx/lcos/wait_all.hpp>

typedef nast_hpx::grid::server::partition_server partition_component;
typedef hpx::components::component<partition_component> partition_server_type;

HPX_REGISTER_COMPONENT_MODULE();

HPX_REGISTER_COMPONENT(partition_server_type, partition_component);

HPX_REGISTER_ACTION(nast_hpx::grid::server::partition_server::do_timestep_action,
    partition_server_do_timestep_action);
HPX_REGISTER_ACTION(nast_hpx::grid::server::partition_server::init_action,
    partition_server_init_action);
    
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
    gy_(cfg.gy)
{    
    data_[U] = partition_data<Real>(cells_x_, cells_y_);
    data_[V] = partition_data<Real>(cells_x_, cells_y_);
    data_[F] = partition_data<Real>(cells_x_, cells_y_);
    data_[G] = partition_data<Real>(cells_x_, cells_y_);
    data_[P] = partition_data<Real>(cells_x_, cells_y_);
    
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
        recv_buffer_top_[P].valid_= true;
        recv_buffer_top_[U].valid_= true;
        recv_buffer_top_[V].valid_= true;
    }
    
    for (uint i = 0; i < data_[P].size_x_ ; i++)
        for (uint j = 0; j < data_[P].size_y_; ++j)            
        {
        //    data_[U](i, j) = 100*(idy_*num_partitions_x_ + idx_)+i+0.01*j;
          /*  data_[V](i, j) = 100*(idy_*num_partitions_x_ + idx_)+i+0.01*j;
            data_[F](i, j) = 100*(idy_*num_partitions_x_ + idx_)+i+0.01*j;
            data_[G](i, j) = 100*(idy_*num_partitions_x_ + idx_)+i+0.01*j;
            data_[P](i, j) = 100*(idy_*num_partitions_x_ + idx_)+i+0.01*j;*/
        }
} 

template<>
void partition_server::send_boundary<LEFT>(std::size_t step, std::vector<hpx::shared_future<void> >& send_futures, std::size_t var)
{
    if (!send_futures.empty())
    {
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
}

template<>
void partition_server::send_boundary<RIGHT>(std::size_t step, std::vector<hpx::shared_future<void> >& send_futures, std::size_t var)
{
    if (!send_futures.empty())
    {
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
}

template<>
void partition_server::send_boundary<BOTTOM>(std::size_t step, std::vector<hpx::shared_future<void> >& send_futures, std::size_t var)
{
    if (!send_futures.empty())
    {
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
}

template<>
void partition_server::send_boundary<TOP>(std::size_t step, std::vector<hpx::shared_future<void> >& send_futures, std::size_t var)
{
    if (!send_futures.empty())
    {
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
}

template<>
void partition_server::send_boundary<BOTTOM_RIGHT>(std::size_t step, std::vector<hpx::shared_future<void> >& send_futures, std::size_t var)
{
    if (!send_futures.empty())
    {
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
}

template<>
void partition_server::send_boundary<TOP_LEFT>(std::size_t step, std::vector<hpx::shared_future<void> >& send_futures, std::size_t var)
{
    //TODO maybe move if to send_?_boundary
    if (!send_futures.empty())
    {
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
}

template<std::size_t var>
void partition_server::send_right_and_top_boundaries(std::size_t step, std::vector<std::vector<hpx::shared_future<void> > >& send_futures)
{
    send_boundary<BOTTOM>(step, send_futures[BOTTOM], var);
    send_boundary<TOP>(step, send_futures[TOP], var);
}

template<std::size_t var>
void partition_server::send_cross_boundaries(std::size_t step, std::vector<std::vector<hpx::shared_future<void> > >& send_futures)
{
    send_boundary<LEFT>(step, send_futures[LEFT], var);
    send_boundary<RIGHT>(step, send_futures[RIGHT], var);
    send_boundary<BOTTOM>(step, send_futures[BOTTOM], var);
    send_boundary<TOP>(step, send_futures[TOP], var);
}

template<std::size_t var>
void partition_server::send_all_boundaries(std::size_t step, std::vector<std::vector<hpx::shared_future<void> > >& send_futures)
{
    send_boundary<LEFT>(step, send_futures[LEFT], var);
    send_boundary<RIGHT>(step, send_futures[RIGHT], var);
    send_boundary<BOTTOM>(step, send_futures[BOTTOM], var);
    send_boundary<TOP>(step, send_futures[TOP], var);
    send_boundary<BOTTOM_RIGHT>(step, send_futures[BOTTOM_RIGHT], var);
    send_boundary<TOP_LEFT>(step, send_futures[TOP_LEFT], var);
}

template<>
hpx::shared_future<void> partition_server::receive_boundary<LEFT>(std::size_t step, std::size_t var)
{
    if (recv_buffer_left_[var].valid_)   
        return hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_left_[var]),
                    boost::ref(data_[var]),
                    step
                )
            );
    else
        return hpx::make_ready_future();
}

template<>
hpx::shared_future<void> partition_server::receive_boundary<RIGHT>(std::size_t step, std::size_t var)
{
    if (recv_buffer_right_[var].valid_)   
        return hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_right_[var]),
                    boost::ref(data_[var]),
                    step
                )
            );
    else
        return hpx::make_ready_future();
}

template<>
hpx::shared_future<void> partition_server::receive_boundary<BOTTOM>(std::size_t step, std::size_t var)
{
    if (recv_buffer_bottom_[var].valid_)   
        return hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_bottom_[var]),
                    boost::ref(data_[var]),
                    step
                )
            );
    else
        return hpx::make_ready_future();
}

template<>
hpx::shared_future<void> partition_server::receive_boundary<TOP>(std::size_t step, std::size_t var)
{
    if (recv_buffer_top_[var].valid_)   
        return hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_top_[var]),
                    boost::ref(data_[var]),
                    step
                )
            );
    else
        return hpx::make_ready_future();
}

template<>
hpx::shared_future<void> partition_server::receive_boundary<BOTTOM_RIGHT>(std::size_t step, std::size_t var)
{
    if (recv_buffer_bottom_right_[var].valid_)   
        return hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_bottom_right_[var]),
                    boost::ref(data_[var]),
                    step
                )
            );
    else
        return hpx::make_ready_future();
}

template<>
hpx::shared_future<void> partition_server::receive_boundary<TOP_LEFT>(std::size_t step, std::size_t var)
{
    if (recv_buffer_top_left_[var].valid_)   
        return hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_top_left_[var]),
                    boost::ref(data_[var]),
                    step
                )
            );
    else
        return hpx::make_ready_future();
}
    
template<std::size_t var>
void partition_server::receive_left_and_bottom_boundaries(std::size_t step, std::vector<hpx::shared_future<void> >& recv_futures)
{
    recv_futures[LEFT] = receive_boundary<LEFT>(step, var);
    recv_futures[BOTTOM] = receive_boundary<BOTTOM>(step, var);
    recv_futures[RIGHT] = hpx::make_ready_future();
    recv_futures[TOP] = hpx::make_ready_future();
    recv_futures[BOTTOM_RIGHT] = hpx::make_ready_future();
    recv_futures[TOP_LEFT] = hpx::make_ready_future();
}  

template<std::size_t var>
void partition_server::receive_cross_boundaries(std::size_t step, std::vector<hpx::shared_future<void> >& recv_futures)
{
    recv_futures[LEFT] = receive_boundary<LEFT>(step, var);
    recv_futures[RIGHT] = receive_boundary<RIGHT>(step, var);
    recv_futures[BOTTOM] = receive_boundary<BOTTOM>(step, var);
    recv_futures[TOP] = receive_boundary<TOP>(step, var);
    recv_futures[BOTTOM_RIGHT] = hpx::make_ready_future();
    recv_futures[TOP_LEFT] = hpx::make_ready_future();
}    

template<std::size_t var>
void partition_server::receive_all_boundaries(std::size_t step, std::vector<hpx::shared_future<void> >& recv_futures)
{
    recv_futures[LEFT] = receive_boundary<LEFT>(step, var);
    recv_futures[RIGHT] = receive_boundary<RIGHT>(step, var);
    recv_futures[BOTTOM] = receive_boundary<BOTTOM>(step, var);
    recv_futures[TOP] = receive_boundary<TOP>(step, var);
    recv_futures[BOTTOM_RIGHT] = receive_boundary<BOTTOM_RIGHT>(step, var);
    recv_futures[TOP_LEFT] = receive_boundary<TOP_LEFT>(step, var);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<LEFT>(std::size_t idx_block,
    std::size_t idy_block, std::vector<hpx::shared_future<void> > const& recv_futures,
    partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if (idx_block == 0)
        return recv_futures[LEFT];
        
    return calc_futures(idx_block - 1, idy_block);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<RIGHT>(std::size_t idx_block,
    std::size_t idy_block, std::vector<hpx::shared_future<void> > const& recv_futures,
    partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if (idx_block == calc_futures.size_x_ - 1)
        return recv_futures[RIGHT];
        
    return calc_futures(idx_block + 1, idy_block);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<BOTTOM>(std::size_t idx_block,
    std::size_t idy_block, std::vector<hpx::shared_future<void> > const& recv_futures,
    partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if (idy_block == 0)
        return recv_futures[BOTTOM];
        
    return calc_futures(idx_block, idy_block - 1);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<TOP>(std::size_t idx_block,
    std::size_t idy_block, std::vector<hpx::shared_future<void> > const& recv_futures,
    partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if (idy_block == calc_futures.size_y_ - 1)
        return recv_futures[TOP];
        
    std::cout << "return calc" << std::endl;
    return calc_futures(idx_block, idy_block + 1);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<BOTTOM_RIGHT>(std::size_t idx_block,
    std::size_t idy_block, std::vector<hpx::shared_future<void> > const& recv_futures,
    partition_data<hpx::shared_future<void> > const& calc_futures)
{
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
    std::size_t idy_block, std::vector<hpx::shared_future<void> > const& recv_futures,
    partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if (idy_block == calc_futures.size_y_ - 1 && idx_block == 0)
        return recv_futures[TOP_LEFT];
    else if (idy_block == calc_futures.size_y_ - 1)
        return recv_futures[TOP];
    else if (idx_block == 0)
        return recv_futures[LEFT];
        
    return calc_futures(idx_block - 1, idy_block + 1);
}
    
void partition_server::do_timestep(Real dt)
{          
    hpx::util::high_resolution_timer t;
           
    partition_data<hpx::shared_future<void> > set_velocity_futures(num_x_blocks_, num_y_blocks_);
    std::vector<std::vector<hpx::shared_future<void> > > send_futures(NUM_DIRECTIONS);
        
    for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += cells_per_y_block_, ++ny_block)
    {
        range_type y_range(y, std::min(y + cells_per_y_block_, cells_y_ - 1));
    
        for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += cells_per_x_block_, ++nx_block)
        {
            range_type x_range(x, std::min(x + cells_per_x_block_, cells_x_ - 1));
                
            std::cout << "block " << nx_block << " " << ny_block << " "
                << x_range.first << " " <<x_range.second <<" " <<y_range.first
                << " " << y_range.second << std::endl;
                
                
            set_velocity_futures(nx_block, ny_block) =
                hpx::async(
                    hpx::util::bind(
                        &stencils<STENCIL_SET_VELOCITY>::call,
                        boost::ref(data_[U]), boost::ref(data_[V]),
                        boost::ref(data_[U]), boost::ref(data_[V]),
                        boost::ref(cell_types_), u_bnd_, v_bnd_,
                        x_range, y_range             
                    )
                );

            if (!is_left_ && nx_block == 0)
                send_futures[LEFT].push_back(set_velocity_futures(nx_block, ny_block));
                
            if (!is_right_ && nx_block == num_x_blocks_ - 1)
                send_futures[RIGHT].push_back(set_velocity_futures(nx_block, ny_block)); 
               
            if (!is_bottom_ && ny_block == 0)
                send_futures[BOTTOM].push_back(set_velocity_futures(nx_block, ny_block));       
        
            if (!is_top_ && ny_block == num_y_blocks_ - 1)
                send_futures[TOP].push_back(set_velocity_futures(nx_block, ny_block));
                
            if (!is_top_ && !is_left_ && nx_block == 0 && ny_block == num_y_blocks_ - 1)
                send_futures[TOP_LEFT].push_back(set_velocity_futures(nx_block, ny_block));
                
            if (!is_bottom_ && !is_right_ && nx_block == num_x_blocks_ - 1 && ny_block == 0)
                send_futures[BOTTOM_RIGHT].push_back(set_velocity_futures(nx_block, ny_block));     
        }
    }
    
    
    std::vector<hpx::shared_future<void> > recv_futures_U(NUM_DIRECTIONS);
    std::vector<hpx::shared_future<void> > recv_futures_V(NUM_DIRECTIONS);
    
    if (rank_ == 3)
    {
        std::cout << "starting wait" << std::endl;
        boost::system_time xt(boost::get_system_time() +
                boost::posix_time::seconds(2));

            boost::thread::sleep(xt);
        std::cout << "end wait" << std::endl;
    }
    send_all_boundaries<U>(0, send_futures);
    send_all_boundaries<V>(0, send_futures);
    
    send_futures.clear();
    
    receive_all_boundaries<U>(0, recv_futures_U);
    receive_all_boundaries<V>(0, recv_futures_V);
    
    if (rank_ == 1)
    {
        std::cout << "waiting for top left" << std::endl;
        get_dependency<TOP>(0, 1, recv_futures_U, set_velocity_futures).wait();
        std::cout << "got it" << std::endl;
    }
    
    
    std::cout << "oooooooooooooooooother work" << std::endl;
    
    for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += cells_per_y_block_, ++ny_block)
    {
        range_type y_range(std::max(y, static_cast<std::size_t>(1 + is_bottom_)), std::min(y + cells_per_y_block_, cells_y_ - 1 - is_top_));
    
        for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += cells_per_x_block_, ++nx_block)
        {
            range_type x_range(std::max(x, static_cast<std::size_t>(1 + is_left_)), std::min(x + cells_per_x_block_, cells_x_ - 1 - is_right_));
            
            std::cout << "block " << nx_block << " " << ny_block << " "
                << x_range.first << " " <<x_range.second <<" " <<y_range.first
                << " " << y_range.second << std::endl;
        }
    }
    
   // hpx::wait_all(recv_futures_U);
   // hpx::wait_all(recv_futures_V);
    hpx::wait_all(set_velocity_futures.data_);

    std::cout << t.elapsed() << std::endl;

    std::cout << "u after\n" << data_[U] << "\n v after\n" << data_[V] << std::endl;
    
    
    for (uint y = cells_y_ - 1; y < cells_y_; --y)
    {
        for (uint x = 0; x < cells_x_; ++x)
            std::cout << cell_types_(x, y).to_ulong() << " ";
        std::cout << "\n";            
    }
    std::cout << std::endl;
    

   
}   
    
}
}
}