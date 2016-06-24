#include <hpx/lcos/gather.hpp>
#include <hpx/lcos/broadcast.hpp>

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

partition_server::partition_server(io::config const& cfg, std::size_t idx, std::size_t idy, std::size_t local_idx, std::size_t local_idy, std::size_t rank)
:   c(cfg),
    cells_x_(c.cells_x_per_partition + 2),
    cells_y_(c.cells_y_per_partition + 2),
    idx_(idx),
    idy_(idy),
    local_idx_(local_idx),
    local_idy_(local_idy),
    rank_(rank),
    is_left_(idx == 0),
    is_right_(idx == c.num_partitions_x - 1),
    is_bottom_(idy == 0),
    is_top_(idy == c.num_partitions_y - 1),
    step_(0),
    next_out_(-1e-12),
    outcount_(0),
    current(1),
    last(0)
{
#ifdef WITH_SOR
    std::cout << "Solver: SOR" << std::endl;
#else
    std::cout << "Solver: blockwise Jacobi" << std::endl;
#endif

#ifdef WITH_FOR_EACH
    std::cout << "Parallelization: hpx::parallel::for_each" << std::endl;
#else
    std::cout << "Parellelization: custom grain size" << std::endl;
#endif

    std::cout << c << std::endl;

    data_[U].resize(cells_x_, cells_y_);
    data_[V].resize(cells_x_, cells_y_);
    data_[F].resize(cells_x_, cells_y_);
    data_[G].resize(cells_x_, cells_y_);
    data_[P].resize(cells_x_, cells_y_);
    rhs_data_.resize(cells_x_, cells_y_);
    cell_type_data_.resize(cells_x_,cells_y_);

    for (std::size_t y = cells_y_ - 1; y < cells_y_; --y)
    {
        for (std::size_t x = 0; x < cells_x_; ++x)
        {
            cell_type_data_(x, y) = std::move(c.flag_grid[(y + local_idy_ * c.cells_y_per_partition) * (c.cells_x_per_partition * c.num_local_partitions_x + 2) + c.cells_x_per_partition * local_idx_ + x]);

            if (x != 0 && x < cells_x_ - 1 && y != 0 && y != cells_y_ - 1)
                if (cell_type_data_(x, y).test(is_fluid))
                    fluid_cells_.emplace_back(x, y);
                else if (cell_type_data_(x, y).test(is_boundary))
                    boundary_cells_.emplace_back(x, y);
                else if (!cell_type_data_(x, y).none())
                    obstacle_cells_.emplace_back(x, y);
        }
            }
   // c.flag_grid.clear();
    
    
   /* for (std::size_t idy_block = 0; idy_block < c.num_y_blocks; ++idy_block)
        for (std::size_t idx_block = 0; idx_block < c.num_x_blocks; ++idx_block)
            std::cout << "BLOCK " << idx_block << " " << idy_block
            << "\nnum_fluid = " << fluid_cells_(idx_block, idy_block).size()
            << "\nnum_boundary = " << boundary_cells_(idx_block, idy_block).size()
            << "\nnum_obstacle = " << obstacle_cells_(idx_block, idy_block).size()
            << std::endl;*/
}

partition_server::partition_server(partition_server && other)
: c(other.c),
  is_left_(other.is_left_),
  is_right_(other.is_right_),
  is_bottom_(other.is_bottom_),
  is_top_(other.is_top_),
  current(other.current),
  last(other.last),
  rhs_data_(std::move(other.rhs_data_)),
  cell_type_data_(other.cell_type_data_),
  fluid_cells_(other.fluid_cells_),
  boundary_cells_(other.boundary_cells_),
  obstacle_cells_(other.obstacle_cells_),
  cells_x_(other.cells_x_),
  cells_y_(other.cells_y_),
  idx_(other.idx_),
  idy_(other.idy_),
  local_idx_(other.local_idx_),
  local_idy_(other.local_idy_),
  rank_(other.rank_),      
  step_(other.step_),
  outcount_(other.outcount_),
  t_(other.t_),
  next_out_(other.next_out_)
{
    for (std::size_t var = 0; var < NUM_VARIABLES; ++var)
        data_[var] = std::move(other.data_[var]);
}


void partition_server::init()
{
    std::vector<hpx::future<hpx::id_type > > parts =
        hpx::find_all_from_basename(partition_basename, c.num_partitions);

    ids_ = hpx::when_all(parts).then(hpx::util::unwrapped2(
                [](std::vector<hpx::id_type>&& ids) -> std::vector<hpx::id_type>
                { return ids;})
            ).get();
    
    if (!is_left_)
    {
        send_buffer_left_.dest_ = ids_[idy_ * c.num_partitions_x + idx_ - 1];
        
        recv_buffer_left_[P].valid_ = true;
        recv_buffer_left_[F].valid_ = true;
        recv_buffer_left_[G].valid_ = true;
        recv_buffer_left_[U].valid_ = true;
        recv_buffer_left_[V].valid_ = true;


        if (!is_top_)
        {           
            send_buffer_top_left_.dest_ = ids_[(idy_ + 1) * c.num_partitions_x + idx_ - 1];

            recv_buffer_top_left_[U].valid_ = true;
            recv_buffer_top_left_[V].valid_ = true;
        }
    }

    if (!is_right_)
    {      
        send_buffer_right_.dest_ = ids_[idy_ * c.num_partitions_x + idx_ + 1];

        recv_buffer_right_[P].valid_ = true;
        recv_buffer_right_[U].valid_ = true;
        recv_buffer_right_[V].valid_ = true;

        if (!is_bottom_)
        {
            send_buffer_bottom_right_.dest_ = ids_[(idy_ - 1) * c.num_partitions_x + idx_ + 1];
            recv_buffer_bottom_right_[U].valid_ = true;
            recv_buffer_bottom_right_[V].valid_ = true;
        }
    }

    if (!is_bottom_)
    {        
        send_buffer_bottom_.dest_ = ids_[(idy_ - 1) * c.num_partitions_x + idx_];

        recv_buffer_bottom_[P].valid_= true;
        recv_buffer_bottom_[F].valid_= true;
        recv_buffer_bottom_[G].valid_= true;
        recv_buffer_bottom_[U].valid_= true;
        recv_buffer_bottom_[V].valid_= true;
    }

    if (!is_top_)
    {
        send_buffer_top_.dest_ = ids_[(idy_ + 1) * c.num_partitions_x + idx_];

        recv_buffer_top_[P].valid_ = true;
        recv_buffer_top_[U].valid_ = true;
        recv_buffer_top_[V].valid_ = true;
    }

    token.reset();
}

template<>
void partition_server::send_boundary<LEFT>(std::size_t step, std::size_t var, hpx::shared_future<void>& send_future)
{
    send_future.then(
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
void partition_server::send_boundary<RIGHT>(std::size_t step, std::size_t var, hpx::shared_future<void>& send_future)
{
    send_future.then(
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
void partition_server::send_boundary<BOTTOM>(std::size_t step, std::size_t var, hpx::shared_future<void>& send_future)
{
    send_future.then(
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
void partition_server::send_boundary<TOP>(std::size_t step, std::size_t var, hpx::shared_future<void>& send_future)
{
    send_future.then(
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
void partition_server::send_boundary<BOTTOM_RIGHT>(std::size_t step, std::size_t var, hpx::shared_future<void>& send_future)
{
    send_future.then(
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
void partition_server::send_boundary<TOP_LEFT>(std::size_t step, std::size_t var, hpx::shared_future<void>& send_future)
{
    send_future.then(
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
void partition_server::send_right_and_top_boundaries(std::size_t step, hpx::shared_future<void>& send_future)
{
    if (!is_right_)
        send_boundary<RIGHT>(step, var, send_future);
        
    if (!is_top_)
        send_boundary<TOP>(step, var, send_future);
}

template<std::size_t var>
void partition_server::send_cross_boundaries(std::size_t step, hpx::shared_future<void>& send_future)
{
    if (!is_left_)
        send_boundary<LEFT>(step, var, send_future);
        
    if (!is_right_)
        send_boundary<RIGHT>(step, var, send_future);
    
    if (!is_bottom_)
        send_boundary<BOTTOM>(step, var, send_future);
    
    if (!is_top_)
        send_boundary<TOP>(step, var, send_future);
}

template<std::size_t var>
void partition_server::send_all_boundaries(std::size_t step, hpx::shared_future<void>& send_future)
{
    if (!is_left_)
        send_boundary<LEFT>(step, var, send_future);
        
    if (!is_right_)
        send_boundary<RIGHT>(step, var, send_future);
    
    if (!is_bottom_)
        send_boundary<BOTTOM>(step, var, send_future);
    
    if (!is_top_)
        send_boundary<TOP>(step, var, send_future);
        
    if (!is_bottom_ && !is_right_)
        send_boundary<BOTTOM_RIGHT>(step, var, send_future);
    
    if (!is_top_ && !is_left_)
        send_boundary<TOP_LEFT>(step, var, send_future);
}


template<>
void partition_server::receive_boundary<LEFT>(std::size_t step, std::size_t var, hpx::shared_future<void>& recv_future)
{   
    recv_future =
        hpx::async(
            hpx::util::bind(
                boost::ref(recv_buffer_left_[var]),
                boost::ref(data_[var]),
                step
            )
    );
}

template<>
void partition_server::receive_boundary<RIGHT>(std::size_t step, std::size_t var, hpx::shared_future<void>& recv_future)
{
    recv_future =
        hpx::async(
            hpx::util::bind(
                boost::ref(recv_buffer_right_[var]),
                boost::ref(data_[var]),
                step
            )
    );
}

template<>
void partition_server::receive_boundary<BOTTOM>(std::size_t step, std::size_t var, hpx::shared_future<void>& recv_future)
{
    recv_future =
        hpx::async(
            hpx::util::bind(
                boost::ref(recv_buffer_bottom_[var]),
                boost::ref(data_[var]),
                step
            )
    );
}

template<>
void partition_server::receive_boundary<TOP>(std::size_t step, std::size_t var, hpx::shared_future<void>& recv_future)
{
    recv_future =
        hpx::async(
            hpx::util::bind(
                boost::ref(recv_buffer_top_[var]),
                boost::ref(data_[var]),
                step
            )
    );
}

template<>
void partition_server::receive_boundary<BOTTOM_RIGHT>(std::size_t step, std::size_t var, hpx::shared_future<void>& recv_future)
{
    recv_future =
        hpx::async(
            hpx::util::bind(
                boost::ref(recv_buffer_bottom_right_[var]),
                boost::ref(data_[var]),
                step
            )
    );
}

template<>
void partition_server::receive_boundary<TOP_LEFT>(std::size_t step, std::size_t var, hpx::shared_future<void>& recv_future)
{
    recv_future =
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
    if (!is_left_)
        receive_boundary<LEFT>(step, var, recv_futures[LEFT]);
    else
        recv_futures[LEFT] = hpx::make_ready_future();
    
    if (!is_bottom_)
        receive_boundary<BOTTOM>(step, var, recv_futures[BOTTOM]);
    else
        recv_futures[BOTTOM] = hpx::make_ready_future(),


    recv_futures[RIGHT] = hpx::make_ready_future();
    recv_futures[TOP] = hpx::make_ready_future();
    recv_futures[BOTTOM_RIGHT] = hpx::make_ready_future();
    recv_futures[TOP_LEFT] = hpx::make_ready_future();
}

template<std::size_t var>
void partition_server::receive_cross_boundaries(std::size_t step, future_vector& recv_futures)
{
    
    if (!is_left_)
        receive_boundary<LEFT>(step, var, recv_futures[LEFT]);
    else
        recv_futures[LEFT] = hpx::make_ready_future();
    
    if (!is_right_)
        receive_boundary<RIGHT>(step, var, recv_futures[RIGHT]);
    else
        recv_futures[RIGHT] = hpx::make_ready_future();
    
    if (!is_bottom_)
        receive_boundary<BOTTOM>(step, var, recv_futures[BOTTOM]);
    else
        recv_futures[BOTTOM] = hpx::make_ready_future();
        
    if (!is_top_)
        receive_boundary<TOP>(step, var, recv_futures[TOP]);
    else
        recv_futures[TOP] = hpx::make_ready_future();

    recv_futures[BOTTOM_RIGHT] = hpx::make_ready_future();
    recv_futures[TOP_LEFT] = hpx::make_ready_future();
}

template<std::size_t var>
void partition_server::receive_all_boundaries(std::size_t step, future_vector& recv_futures)
{
    if (!is_left_)
        receive_boundary<LEFT>(step, var, recv_futures[LEFT]);
    else
        recv_futures[LEFT] = hpx::make_ready_future();
    
    if (!is_right_)
        receive_boundary<RIGHT>(step, var, recv_futures[RIGHT]);
    else
        recv_futures[RIGHT] = hpx::make_ready_future();
    
    if (!is_bottom_)
        receive_boundary<BOTTOM>(step, var, recv_futures[BOTTOM]);
    else
        recv_futures[BOTTOM] = hpx::make_ready_future();
        
    if (!is_top_)
        receive_boundary<TOP>(step, var, recv_futures[TOP]);
    else
        recv_futures[TOP] = hpx::make_ready_future();

    if (!is_bottom_ && !is_right_)
        receive_boundary<BOTTOM_RIGHT>(step, var, recv_futures[BOTTOM_RIGHT]);
    else
        recv_futures[BOTTOM_RIGHT] = hpx::make_ready_future();
        
    if (!is_top_ && !is_left_)
        receive_boundary<TOP_LEFT>(step, var, recv_futures[TOP_LEFT]);
    else
        recv_futures[TOP_LEFT] = hpx::make_ready_future();
}

void partition_server::wait_all_boundaries(future_vector& recv_futures)
{
    if (!is_left_)
        recv_futures[LEFT].wait();

    if (!is_right_)
        recv_futures[RIGHT].wait();

    if (!is_bottom_)
    {
        recv_futures[BOTTOM].wait();

        if (!is_right_)
            recv_futures[BOTTOM_RIGHT].wait();
    }

    if (!is_top_)
    {
        recv_futures[TOP].wait();

        if (!is_left_)
            recv_futures[TOP_LEFT].wait();
    }
}

std::pair<Real, Real> partition_server::do_timestep(Real dt)
{
    hpx::util::high_resolution_timer t;
    
    hpx::shared_future<void> set_velocity_future1 =
        hpx::async(
            hpx::util::bind(
                &stencils<STENCIL_SET_VELOCITY_BOUNDARY>::call,
                boost::ref(data_[U]), boost::ref(data_[V]),
                boost::ref(data_[U]), boost::ref(data_[V]),
                boost::ref(cell_type_data_), boost::ref(boundary_cells_),
                c.u_bnd, c.v_bnd, c.bnd_type
            )
        );

    hpx::shared_future<void> set_velocity_future2 =
        hpx::async(
            hpx::util::bind(
                &stencils<STENCIL_SET_VELOCITY_OBSTACLE>::call,
                boost::ref(data_[U]), boost::ref(data_[V]),
                boost::ref(data_[U]), boost::ref(data_[V]),
                boost::ref(cell_type_data_), boost::ref(obstacle_cells_)
            )
        );

    hpx::shared_future<void> set_velocity_future3 =
        hpx::async(
            hpx::util::bind(
                &stencils<STENCIL_SET_VELOCITY_FLUID>::call,
                boost::ref(data_[U]), boost::ref(data_[V]),
                boost::ref(data_[U]), boost::ref(data_[V]),
                boost::ref(cell_type_data_), boost::ref(fluid_cells_)
            )
        );

    hpx::shared_future<void> set_velocity_future =
        hpx::when_all(set_velocity_future1, set_velocity_future2, set_velocity_future3).then(
            [](hpx::lcos::future<hpx::util::tuple<hpx::lcos::shared_future<void>, hpx::lcos::shared_future<void>, hpx::lcos::shared_future<void> > >)
            {return;});

    send_all_boundaries<U>(step_, set_velocity_future);
    
    send_all_boundaries<V>(step_, set_velocity_future);

    future_vector recv_futures_U(NUM_DIRECTIONS);
    future_vector recv_futures_V(NUM_DIRECTIONS);       

    receive_all_boundaries<U>(step_, recv_futures_U);
    receive_all_boundaries<V>(step_, recv_futures_V);
        
    if (c.vtk && next_out_ < t_)
    {
        next_out_ += c.delta_vec;

        std::cout << "Output in step " << step_ << " " << c.eps << " " << c.iter_max << token.was_cancelled() << std::endl;

        set_velocity_future.then(
            hpx::launch::async,
            hpx::util::bind(
                &io::writer::write_vtk,
                boost::ref(data_[P]), boost::ref(data_[U]), boost::ref(data_[V]), boost::ref(cell_type_data_),
                c.num_partitions_x, c.num_partitions_y, c.i_max, c.j_max, c.dx, c.dx, outcount_++,
                c.rank
            )
        ).wait();
    }

 

    hpx::shared_future<void> compute_fg_future1 =
        hpx::dataflow(
            hpx::util::unwrapped(
                hpx::util::bind(
                    &stencils<STENCIL_COMPUTE_FG_BOUNDARY_AND_OBSTACLE>::call,
                    boost::ref(data_[F]), boost::ref(data_[G]),
                    boost::ref(data_[U]), boost::ref(data_[V]),
                    boost::ref(cell_type_data_),
                    boost::ref(boundary_cells_),
                    boost::ref(obstacle_cells_)
                )
            )
            , set_velocity_future
            , recv_futures_U[LEFT]
            , recv_futures_V[LEFT]
            , recv_futures_U[RIGHT]
            , recv_futures_V[RIGHT]
            , recv_futures_U[BOTTOM]
            , recv_futures_V[BOTTOM]
            , recv_futures_U[TOP]
            , recv_futures_V[TOP]
            , recv_futures_U[TOP_LEFT]
            , recv_futures_V[BOTTOM_RIGHT]
        );
    
    hpx::shared_future<void> compute_fg_future2 =
        hpx::dataflow(
            hpx::util::unwrapped(
                hpx::util::bind(
                    &stencils<STENCIL_COMPUTE_FG_FLUID>::call,
                    boost::ref(data_[F]), boost::ref(data_[G]),
                    boost::ref(data_[U]), boost::ref(data_[V]),
                    boost::ref(cell_type_data_),
                    boost::ref(fluid_cells_),
                    c.re, c.gx, c.gy, c.beta, c.dx, c.dy, dt, c.alpha
                )
            )
            , set_velocity_future
            , recv_futures_U[LEFT]
            , recv_futures_V[LEFT]
            , recv_futures_U[RIGHT]
            , recv_futures_V[RIGHT]
            , recv_futures_U[BOTTOM]
            , recv_futures_V[BOTTOM]
            , recv_futures_U[TOP]
            , recv_futures_V[TOP]
            , recv_futures_U[TOP_LEFT]
            , recv_futures_V[BOTTOM_RIGHT]
        );

    hpx::shared_future<void> compute_fg_future =
        hpx::when_all(compute_fg_future1, compute_fg_future2).then(
            [](hpx::lcos::future<hpx::util::tuple<hpx::lcos::shared_future<void>, hpx::lcos::shared_future<void>> >)
            {return;});
    
    send_right_and_top_boundaries<F>(step_, compute_fg_future);
    send_right_and_top_boundaries<G>(step_, compute_fg_future);

    future_vector recv_futures_F(NUM_DIRECTIONS);
    future_vector recv_futures_G(NUM_DIRECTIONS);
    
    receive_left_and_bottom_boundaries<F>(step_, recv_futures_F);
    receive_left_and_bottom_boundaries<G>(step_, recv_futures_G);

    hpx::shared_future<void> compute_rhs_future =
        hpx::dataflow(
            hpx::util::unwrapped(
                hpx::util::bind(
                    &stencils<STENCIL_COMPUTE_RHS>::call,
                    boost::ref(rhs_data_), boost::ref(data_[F]),
                    boost::ref(data_[G]), boost::ref(cell_type_data_),
                    boost::ref(fluid_cells_),
                    c.dx, c.dy, dt
                )
            )
            , compute_fg_future
            , recv_futures_F[LEFT]
            , recv_futures_G[BOTTOM]
        );
   
    hpx::util::high_resolution_timer t1;

    Real rres;
    token.reset();
    
    std::array<hpx::shared_future<void>, 2> sor_cycle_futures;
    std::array<future_vector, 2> recv_futures_P;
    
    for (auto& vec : recv_futures_P)
        vec.resize(NUM_DIRECTIONS);
        
    for (std::size_t dir = 0; dir < NUM_DIRECTIONS; ++dir)
        recv_futures_P[last][dir] = hpx::make_ready_future();
    
    sor_cycle_futures[last] = hpx::make_ready_future();
    
    hpx::shared_future<Real> compute_res_future = hpx::make_ready_future(0.);

    for (std::size_t iter = 0; iter < c.iter_max; ++iter)
    {
#ifdef WITH_SOR
       
        hpx::shared_future<void> set_p_future =
            hpx::dataflow(
                hpx::util::unwrapped(
                    hpx::util::bind(
                        &stencils<STENCIL_SET_P>::call,
                        boost::ref(data_[P]), boost::ref(cell_type_data_),
                        boost::ref(boundary_cells_),
                        boost::ref(obstacle_cells_),
                        token
                    )
                )
                , sor_cycle_futures[last]
                , recv_futures_P[last][LEFT]
                , recv_futures_P[last][RIGHT]
                , recv_futures_P[last][BOTTOM]
                , recv_futures_P[last][TOP]
                , compute_res_future.then([](hpx::shared_future<Real> a) {return;})
            );
 

        receive_cross_boundaries<P>(iter, recv_futures_P[current]);

        sor_cycle_futures[current] =
            hpx::dataflow(
                hpx::util::unwrapped(
                    hpx::util::bind(
                        &stencils<STENCIL_SOR>::call,
                        boost::ref(data_[P]), boost::ref(rhs_data_),
                        boost::ref(fluid_cells_),
                        c.part1, c.part2, c.dx_sq, c.dy_sq,
                        token, iter
                    )
                )
                , set_p_future
                , compute_rhs_future
                , recv_futures_P[current][LEFT]
                , recv_futures_P[current][BOTTOM]
                , recv_futures_P[last][RIGHT]
                , recv_futures_P[last][TOP]
            );    

        send_cross_boundaries<P>(iter, sor_cycle_futures[current]);

        compute_res_future =
            hpx::dataflow(
                hpx::util::unwrapped(
                    hpx::util::bind(
                        &stencils<STENCIL_COMPUTE_RESIDUAL>::call,
                        boost::ref(data_[P]), boost::ref(rhs_data_),
                        boost::ref(fluid_cells_),
                        c.over_dx_sq, c.over_dy_sq,
                        token
                    )
                )
                , sor_cycle_futures[current]
                , recv_futures_P[current][RIGHT]
                , recv_futures_P[current][TOP]
            );

#else

        hpx::shared_future<void> set_p_future =
            hpx::dataflow(
                hpx::util::unwrapped(
                    hpx::util::bind(
                        &stencils<STENCIL_SET_P>::call,
                        boost::ref(data_[P]), boost::ref(cell_type_data_),
                        boost::ref(boundary_cells_), boost::ref(obstacle_cells_),
                        token
                    )
                )
                , sor_cycle_futures[last]
                , recv_futures_P[last][LEFT]
                , recv_futures_P[last][RIGHT]
                , recv_futures_P[last][BOTTOM]
                , recv_futures_P[last][TOP]
                , compute_res_future.then([](hpx::shared_future<Real> a) {return;})
            );


        sor_cycle_futures[current] =
            hpx::dataflow(
                hpx::util::unwrapped(
                    hpx::util::bind(
                        &stencils<STENCIL_JACOBI>::call,
                        boost::ref(data_[P]), boost::ref(rhs_data_),
                        boost::ref(fluid_cells_), c.dx_sq, c.dy_sq,
                        token, iter
                    )
                )
                , set_p_future
                , compute_rhs_future
            );

        send_cross_boundaries<P>(iter, sor_cycle_futures[current]);
        receive_cross_boundaries<P>(iter, recv_futures_P[current]);

        compute_res_future =
            hpx::dataflow(
                hpx::util::unwrapped(
                    hpx::util::bind(
                        &stencils<STENCIL_COMPUTE_RESIDUAL>::call,
                        boost::ref(data_[P]), boost::ref(rhs_data_),
                        boost::ref(fluid_cells_), c.over_dx_sq, c.over_dy_sq,
                       token
                    )
                )
                , sor_cycle_futures[current]
                , recv_futures_P[current][LEFT]
                , recv_futures_P[current][RIGHT]
                , recv_futures_P[current][BOTTOM]
                , recv_futures_P[current][TOP]
            );

#endif
        std::swap(current, last);

        hpx::future<Real> local_residual =
            compute_res_future.then([](hpx::shared_future<Real> a) -> Real {return a.get();});

        if (idx_ == 0 && idy_ == 0)
        {          
            hpx::future<std::vector<Real> > partial_residuals =
                hpx::lcos::gather_here(residual_basename,
                                        std::move(local_residual),
                                        c.num_partitions, step_ * c.iter_max + iter, 0);

            hpx::future<Real> residual =
                partial_residuals.then(
                    hpx::util::unwrapped(
                        [this](std::vector<Real> local_residuals)
                            -> Real
                        {
                            Real residual = 0;

                            for (std::size_t i = 0; i < local_residuals.size(); ++i)
                                residual += local_residuals[i];

                            return std::sqrt(residual/c.num_fluid_cells);
                        }
                    )
                );

            residual.then(
                hpx::util::unwrapped(
                    [dt, iter, step = step_, t = t_, this](Real residual)
                    {
                        if ((residual < c.eps || iter == c.iter_max - 1)
                            && !token.was_cancelled())
                        {
                            std::cout << "step = " << step
                                << ", t = " << t
                                << ", dt = " << dt
                                << ", iter = "<< iter
                                << ", residual = " << residual
                                << std::endl;

                            hpx::lcos::broadcast_apply<cancel_action>(ids_);
                        }
                    }
                )
            );
        }
        // if not root locality, send residual to root locality
        else
            hpx::lcos::gather_there(residual_basename, std::move(local_residual),
                                        step_ * c.iter_max + iter, 0, rank_);
    }

    hpx::shared_future<std::pair<Real, Real> > local_max_uv =
        hpx::dataflow(
            hpx::util::unwrapped(
                hpx::util::bind(
                    &stencils<STENCIL_UPDATE_VELOCITY>::call,
                    boost::ref(data_[U]), boost::ref(data_[V]), boost::ref(data_[F]),
                    boost::ref(data_[G]), boost::ref(data_[P]),
                    boost::ref(cell_type_data_),
                    boost::ref(fluid_cells_),
                    dt, c.over_dx, c.over_dy
                )
            )
            , recv_futures_P[last][LEFT]
            , recv_futures_P[last][RIGHT]
            , recv_futures_P[last][BOTTOM]
            , recv_futures_P[last][TOP]
            , sor_cycle_futures[last]
        );

    t_ += dt;
    ++step_;
    
    return local_max_uv.get();
}


}
}
}
