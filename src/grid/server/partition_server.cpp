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

partition_server::partition_server(io::config&& cfg)
:   c(cfg),
    cells_x_(c.cells_x_per_block * c.num_x_blocks + 2),
    cells_y_(c.cells_y_per_block * c.num_y_blocks + 2),
    is_left_(c.idx == 0),
    is_right_(c.idx == c.num_partitions_x - 1),
    is_bottom_(c.idy == 0),
    is_top_(c.idy == c.num_partitions_y - 1),
    step_(0),
    next_out_(-1e-12),
    outcount_(0)
{
    std::cout << c << std::endl;

    data_[U].resize(cells_x_, cells_y_);
    data_[V].resize(cells_x_, cells_y_);
    data_[F].resize(cells_x_, cells_y_);
    data_[G].resize(cells_x_, cells_y_);
    data_[P].resize(cells_x_, cells_y_);
    rhs_data_.resize(cells_x_, cells_y_);
    cell_type_data_.resize(cells_x_,cells_y_);

   for (std::size_t i = 0; i != cell_type_data_.size_; ++i)
        cell_type_data_[i] = std::move(c.flag_grid[i]);

    c.flag_grid.clear();
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
        send_buffer_left_.dest_ = ids_[c.idy * c.num_partitions_x + c.idx - 1];

        recv_buffer_left_[P].valid_ = true;
        recv_buffer_left_[F].valid_ = true;
        recv_buffer_left_[G].valid_ = true;
        recv_buffer_left_[U].valid_ = true;
        recv_buffer_left_[V].valid_ = true;


        if (!is_top_)
        {
            send_buffer_top_left_.dest_ = ids_[(c.idy + 1) * c.num_partitions_x + c.idx - 1];

            recv_buffer_top_left_[U].valid_ = true;
            recv_buffer_top_left_[V].valid_ = true;
        }
    }

    if (!is_right_)
    {
        send_buffer_right_.dest_ = ids_[c.idy * c.num_partitions_x + c.idx + 1];

        recv_buffer_right_[P].valid_ = true;
        recv_buffer_right_[U].valid_ = true;
        recv_buffer_right_[V].valid_ = true;

        if (!is_bottom_)
        {
            send_buffer_bottom_right_.dest_ = ids_[(c.idy - 1) * c.num_partitions_x + c.idx + 1];
            recv_buffer_bottom_right_[U].valid_ = true;
            recv_buffer_bottom_right_[V].valid_ = true;
        }
    }

    if (!is_bottom_)
    {
        send_buffer_bottom_.dest_ = ids_[(c.idy - 1) * c.num_partitions_x + c.idx];

        recv_buffer_bottom_[P].valid_= true;
        recv_buffer_bottom_[F].valid_= true;
        recv_buffer_bottom_[G].valid_= true;
        recv_buffer_bottom_[U].valid_= true;
        recv_buffer_bottom_[V].valid_= true;
    }

    if (!is_top_)
    {
        send_buffer_top_.dest_ = ids_[(c.idy + 1) * c.num_partitions_x + c.idx];

        recv_buffer_top_[P].valid_ = true;
        recv_buffer_top_[U].valid_ = true;
        recv_buffer_top_[V].valid_ = true;
    }

    for (uint i = 0; i < data_[P].size_x_ ; i++)
        for (uint j = 0; j < data_[P].size_y_; ++j)
        {
          //  data_[U](i, j) = 100*(c.idy*c.num_partitions_x + c.idx)+i+0.01*j;
          /*  data_[V](i, j) = 100*(c.idy*c.num_partitions_x + c.idx)+i+0.01*j;
            data_[F](i, j) = 100*(c.idy*c.num_partitions_x + c.idx)+i+0.01*j;
            data_[G](i, j) = 100*(c.idy*c.num_partitions_x + c.idx)+i+0.01*j;
            data_[P](i, j) = 100*(c.idy*c.num_partitions_x + c.idx)+i+0.01*j;*/
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

    partition_data<hpx::shared_future<void> > set_velocity_futures(c.num_x_blocks, c.num_y_blocks);

    future_grid send_futures_U(NUM_DIRECTIONS);
    future_grid send_futures_V(NUM_DIRECTIONS);

   // std::cout << "compute uv" << std::endl;
  //  std::cout << "before U\n" << data_[U] << std::endl;
   // std::cout << "before V\n" << data_[V] << std::endl;

    for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
    {
        range_type y_range(y, std::min(y + c.cells_y_per_block, cells_y_ - 1));

        for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
        {
            range_type x_range(x, std::min(x + c.cells_x_per_block, cells_x_ - 1));

           // std::cout << "block " << nx_block << " " << ny_block << " "
             //   << x_range.first << " " <<x_range.second <<" " <<y_range.first
               // << " " << y_range.second << std::endl;

            hpx::shared_future<void> calc_future =
                hpx::async(
                    hpx::util::bind(
                        &stencils<STENCIL_SET_VELOCITY>::call,
                        boost::ref(data_[U]), boost::ref(data_[V]),
                        boost::ref(data_[U]), boost::ref(data_[V]),
                        boost::ref(cell_type_data_), c.u_bnd, c.v_bnd,
                        x_range, y_range
                    )
                );

            set_velocity_futures(nx_block, ny_block) = calc_future;

            if (!is_left_ && nx_block == 0)
            {
                send_futures_U[LEFT].push_back(calc_future);
                send_futures_V[LEFT].push_back(calc_future);
            }

            if (!is_right_ && nx_block == c.num_x_blocks - 1)
            {
                send_futures_U[RIGHT].push_back(calc_future);
                send_futures_V[RIGHT].push_back(calc_future);
            }

            if (!is_bottom_ && ny_block == 0)
            {
                send_futures_U[BOTTOM].push_back(calc_future);
                send_futures_V[BOTTOM].push_back(calc_future);
            }

            if (!is_top_ && ny_block == c.num_y_blocks - 1)
            {
                send_futures_U[TOP].push_back(calc_future);
                send_futures_V[TOP].push_back(calc_future);
            }

            if (!is_bottom_ && !is_right_ && nx_block == c.num_x_blocks - 1 && ny_block == 0)
            {
                send_futures_U[BOTTOM_RIGHT].push_back(calc_future);
                send_futures_V[BOTTOM_RIGHT].push_back(calc_future);
            }

            if (!is_top_ && !is_left_ && nx_block == 0 && ny_block == c.num_y_blocks - 1)
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

    if (c.vtk && next_out_ < t_)
    {
        std::cout << "output in step " << step_ << std::endl;

        next_out_ += c.delta_vec;

        hpx::when_all(set_velocity_futures.data_).then(
            hpx::launch::async,
            hpx::util::bind(
                &io::writer::write_vtk,
                boost::ref(data_[P]), boost::ref(data_[U]), boost::ref(data_[V]),
                c.num_partitions_x, c.num_partitions_y, c.i_max, c.j_max, c.dx, c.dx, outcount_++,
                c.rank
            )
        ).wait();
    }

    partition_data<hpx::shared_future<void> > compute_fg_futures(c.num_x_blocks, c.num_y_blocks);
  //  std::cout << "compute fg" << std::endl;

    future_vector dependencies;
    future_grid send_futures_F(NUM_DIRECTIONS);
    future_grid send_futures_G(NUM_DIRECTIONS);

    for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
    {
        range_type y_range(y, std::min(y + c.cells_y_per_block, cells_y_ - 1));

        for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
        {
            range_type x_range(x, std::min(x + c.cells_x_per_block, cells_x_ - 1));

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
                        boost::ref(cell_type_data_), c.re, c.gx, c.gy, c.beta, c.dx,
                        c.dy, dt, c.alpha, x_range, y_range
                    )
                );

            compute_fg_futures(nx_block, ny_block) = calc_future;

            if (!is_right_ && nx_block == c.num_x_blocks - 1)
            {
                send_futures_F[RIGHT].push_back(calc_future);
                send_futures_G[RIGHT].push_back(calc_future);
            }

            if (!is_top_ && ny_block == c.num_y_blocks - 1)
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


    partition_data<hpx::shared_future<void> > compute_rhs_futures(c.num_x_blocks, c.num_y_blocks);
  //  std::cout << "compute rhs" << std::endl;

    for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
    {
        range_type y_range(y, std::min(y + c.cells_y_per_block, cells_y_ - 1));

        for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
        {
            range_type x_range(x, std::min(x + c.cells_x_per_block, cells_x_ - 1));

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
                        boost::ref(data_[G]), boost::ref(cell_type_data_),
                        c.dx, c.dy, dt, x_range, y_range
                    )
                );
        }
    }

   /*  hpx::wait_all(compute_rhs_futures.data_);

    std::cout << "RHS\n" << rhs_data_ << std::endl;*/

    //TODO make calc_futures members
    partition_data<hpx::shared_future<void> > set_p_futures(c.num_x_blocks, c.num_y_blocks);
    std::array<partition_data<hpx::shared_future<void> >, 2> sor_cycle_futures;

    std::array<future_vector, 2> recv_futures_P;

    hpx::util::high_resolution_timer t1;

    std::size_t current = 1;
    std::size_t last = 0;

    sor_cycle_futures[current].resize(c.num_x_blocks, c.num_y_blocks);
    sor_cycle_futures[last].resize(c.num_x_blocks, c.num_y_blocks);

    recv_futures_P[current].resize(NUM_DIRECTIONS);
    recv_futures_P[last].resize(NUM_DIRECTIONS);

    for (auto& a : sor_cycle_futures[last].data_)
        a = hpx::make_ready_future();

    for (auto& a : recv_futures_P[last])
        a = hpx::make_ready_future();


    //std::cout << "start sor"<<std::endl;

    Real rres;
    std::size_t iter = 0;
    for (iter = 0; iter < c.iter_max; ++iter)
    {
        std::vector<hpx::future<Real> > local_residuals;
        future_grid send_futures_P(NUM_DIRECTIONS);

       // hpx::wait_all(sor_cycle_futures[last].data_);
        //std::cout << "Before set P\n" << data_[P] << std::endl;


        receive_cross_boundaries<P>(iter, recv_futures_P[current]);

        for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
        {
            range_type y_range(y, std::min(y + c.cells_y_per_block, cells_y_ - 1));

            for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
            {
                range_type x_range(x, std::min(x + c.cells_x_per_block, cells_x_ - 1));

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
                            boost::ref(data_[P]), boost::ref(cell_type_data_),
                            nx_block, ny_block,
                            x_range, y_range
                        )
                    );
            }
        }

      //  hpx::wait_all(set_p_futures.data_);

      //  std::cout << "After set P\n" << data_[P] << std::endl;

        for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
        {
            range_type y_range(y, std::min(y + c.cells_y_per_block, cells_y_ - 1));

            for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
            {
                range_type x_range(x, std::min(x + c.cells_x_per_block, cells_x_ - 1));

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
                            boost::ref(cell_type_data_), c.part1, c.part2, c.dx_sq, c.dy_sq,
                            nx_block, ny_block, iter,
                            x_range, y_range
                        )
                    );

                sor_cycle_futures[current](nx_block, ny_block) = calc_future;

                if (!is_left_ && nx_block == 0)
                    send_futures_P[LEFT].push_back(calc_future);

                if (!is_right_ && nx_block == c.num_x_blocks - 1)
                    send_futures_P[RIGHT].push_back(calc_future);

                if (!is_bottom_ && ny_block == 0)
                    send_futures_P[BOTTOM].push_back(calc_future);

                if (!is_top_ && ny_block == c.num_y_blocks - 1)
                    send_futures_P[TOP].push_back(calc_future);
            }
        }

       // hpx::wait_all(sor_cycle_futures[current].data_);
      //  std::cout << "after sor P\n" << data_[P] << std::endl;

        send_cross_boundaries<P>(iter, send_futures_P);

        for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
        {
            range_type y_range(y, std::min(y + c.cells_y_per_block, cells_y_ - 1));

            for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
            {
                range_type x_range(x, std::min(x + c.cells_x_per_block, cells_x_ - 1));

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
                            boost::ref(cell_type_data_), c.over_dx_sq, c.over_dy_sq,
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
                        [num_fluid_cells = c.num_fluid_cells](std::vector<Real> residuals)
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
                                        c.num_partitions, step_ * c.iter_max + iter);

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
                                        step_ * c.iter_max + iter).wait();
    }



  //  hpx::wait_all(sor_cycle_futures[last].data_);
   //  std::cout << "before2 U\n" << data_[U] << std::endl;
   // std::cout << "before2 V\n" << data_[V] << std::endl;
    Real elapsed = t1.elapsed();

    std::vector<hpx::future<std::pair<Real, Real> > > local_max_uvs;

    for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
    {
        range_type y_range(y, std::min(y + c.cells_y_per_block, cells_y_ - 1));

        for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
        {
            range_type x_range(x, std::min(x + c.cells_x_per_block, cells_x_ - 1));

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
                        boost::ref(cell_type_data_), dt, c.over_dx, c.over_dy,
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
        std::cout << "step " << step_ << " t " << t_ << " dt " << dt << " iter " << iter << " residual " << rres << " elapsed " << elapsed << " average " << elapsed / c.iter_max << std::endl;

    t_ += dt;
    ++step_;

    return m;
}

}
}
}
