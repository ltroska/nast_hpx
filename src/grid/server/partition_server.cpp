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
    outcount_(0),
    send_futures_U(NUM_DIRECTIONS),
    send_futures_V(NUM_DIRECTIONS),
    recv_futures_U(NUM_DIRECTIONS),
    recv_futures_V(NUM_DIRECTIONS),
    send_futures_F(NUM_DIRECTIONS),
    send_futures_G(NUM_DIRECTIONS),
    recv_futures_F(NUM_DIRECTIONS),
    recv_futures_G(NUM_DIRECTIONS),
    current(1),
    last(0),
    send_futures_P(NUM_DIRECTIONS)
{
    std::cout << c << std::endl;

#ifdef WITH_SOR
    std::cout << "Solver: SOR" << std::endl;
#else
    std::cout << "Solver: blockwise Jacobi" << std::endl;
#endif

    data_[U].resize(cells_x_, cells_y_);
    data_[V].resize(cells_x_, cells_y_);
    data_[F].resize(cells_x_, cells_y_);
    data_[G].resize(cells_x_, cells_y_);
    data_[P].resize(cells_x_, cells_y_);
    rhs_data_.resize(cells_x_, cells_y_);
    cell_type_data_.resize(cells_x_,cells_y_);

    fluid_cells_.resize(c.num_x_blocks, c.num_y_blocks);
    boundary_cells_.resize(c.num_x_blocks, c.num_y_blocks);
    obstacle_cells_.resize(c.num_x_blocks, c.num_y_blocks);

    for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
    {
        range_type y_range(y, std::min(y + c.cells_y_per_block, cells_y_ - 1));

        for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
        {
            range_type x_range(x, std::min(x + c.cells_x_per_block, cells_x_ - 1));

            for (std::size_t j = y_range.first; j < y_range.second; ++j)
                for (std::size_t i = x_range.first; i < x_range.second; ++i)
                {
                    cell_type_data_(i, j) = std::move(c.flag_grid[j * cells_x_ + i]);

                    if (cell_type_data_(i, j).test(is_fluid))
                        fluid_cells_(nx_block, ny_block).emplace_back(i, j);
                    else if (cell_type_data_(i, j).test(is_boundary))
                        boundary_cells_(nx_block, ny_block).emplace_back(i, j);
                    else if (!cell_type_data_(i, j).none())
                        obstacle_cells_(nx_block, ny_block).emplace_back(i, j);
                }

        }
    }

    c.flag_grid.clear();

    set_velocity_futures.resize(c.num_x_blocks, c.num_y_blocks);
    init_future_grid(recv_futures_U);
    init_future_grid(recv_futures_V);

    compute_fg_futures.resize(c.num_x_blocks, c.num_y_blocks);
    init_future_grid(recv_futures_F);
    init_future_grid(recv_futures_G);

    compute_rhs_futures.resize(c.num_x_blocks, c.num_y_blocks);
    compute_res_futures.resize(c.num_x_blocks, c.num_y_blocks);

    for (auto& a : compute_res_futures)
        a = hpx::make_ready_future(0.);

    set_p_futures.resize(c.num_x_blocks, c.num_y_blocks);
    recv_futures_P[current].resize(NUM_DIRECTIONS);
    recv_futures_P[last].resize(NUM_DIRECTIONS);
    init_future_grid(recv_futures_P[current]);
    init_future_grid(recv_futures_P[last]);

    for (auto& a : recv_futures_P[last])
        for (auto& b : a)
            b = hpx::make_ready_future();

    sor_cycle_futures[current].resize(c.num_x_blocks, c.num_y_blocks);
    sor_cycle_futures[last].resize(c.num_x_blocks, c.num_y_blocks);

    for (auto& a : sor_cycle_futures[last])
        a = hpx::make_ready_future();

   /* for (std::size_t idy_block = 0; idy_block < c.num_y_blocks; ++idy_block)
        for (std::size_t idx_block = 0; idx_block < c.num_x_blocks; ++idx_block)
            std::cout << "BLOCK " << idx_block << " " << idy_block
            << "\nnum_fluid = " << fluid_cells_(idx_block, idy_block).size()
            << "\nnum_boundary = " << boundary_cells_(idx_block, idy_block).size()
            << "\nnum_obstacle = " << obstacle_cells_(idx_block, idy_block).size()
            << std::endl;*/
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


}

template<>
void partition_server::send_boundary<LEFT>(std::size_t step, std::size_t var, future_vector& send_futures)
{
    if (!send_futures.empty())
        for (std::size_t i = 0; i < send_futures.size(); ++i)
            send_futures[i].then(
                hpx::launch::async,
                hpx::util::bind(
                    boost::ref(send_buffer_left_),
                    boost::ref(data_[var]),
                    step * c.num_y_blocks + i,
                    var,
                    i * c.cells_y_per_block,
                    c.cells_y_per_block
                )
            );
}

template<>
void partition_server::send_boundary<RIGHT>(std::size_t step, std::size_t var, future_vector& send_futures)
{
    if (!send_futures.empty())
        for (std::size_t i = 0; i < send_futures.size(); ++i)
            send_futures[i].then(
                hpx::launch::async,
                hpx::util::bind(
                    boost::ref(send_buffer_right_),
                    boost::ref(data_[var]),
                    step * c.num_y_blocks + i,
                    var,
                    i * c.cells_y_per_block,
                    c.cells_y_per_block
                )
            );
}

template<>
void partition_server::send_boundary<BOTTOM>(std::size_t step, std::size_t var, future_vector& send_futures)
{
    if (!send_futures.empty())
        for (std::size_t i = 0; i < send_futures.size(); ++i)
            send_futures[i].then(
                hpx::launch::async,
                hpx::util::bind(
                    boost::ref(send_buffer_bottom_),
                    boost::ref(data_[var]),
                    step * c.num_x_blocks + i,
                    var,
                    i * c.cells_x_per_block,
                    c.cells_x_per_block
                )
            );
}

template<>
void partition_server::send_boundary<TOP>(std::size_t step, std::size_t var, future_vector& send_futures)
{
    if (!send_futures.empty())
        for (std::size_t i = 0; i < send_futures.size(); ++i)
            send_futures[i].then(
                hpx::launch::async,
                hpx::util::bind(
                    boost::ref(send_buffer_top_),
                    boost::ref(data_[var]),
                    step * c.num_x_blocks + i,
                    var,
                    i * c.cells_x_per_block,
                    c.cells_x_per_block
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
        for (std::size_t i = 0; i < recv_futures.size(); ++i)
            recv_futures[i] =
                hpx::async(
                    hpx::util::bind(
                        boost::ref(recv_buffer_left_[var]),
                        boost::ref(data_[var]),
                        step * c.num_y_blocks + i,
                        i * c.cells_y_per_block
                    )
            );
}

template<>
void partition_server::receive_boundary<RIGHT>(std::size_t step, std::size_t var, future_vector& recv_futures)
{
    if (recv_buffer_right_[var].valid_)
        for (std::size_t i = 0; i < recv_futures.size(); ++i)
            recv_futures[i] =
                hpx::async(
                    hpx::util::bind(
                        boost::ref(recv_buffer_right_[var]),
                        boost::ref(data_[var]),
                        step * c.num_y_blocks + i,
                        i * c.cells_y_per_block
                    )
            );
}

template<>
void partition_server::receive_boundary<BOTTOM>(std::size_t step, std::size_t var, future_vector& recv_futures)
{
    if (recv_buffer_bottom_[var].valid_)
        for (std::size_t i = 0; i < recv_futures.size(); ++i)
            recv_futures[i] =
                hpx::async(
                    hpx::util::bind(
                        boost::ref(recv_buffer_bottom_[var]),
                        boost::ref(data_[var]),
                        step * c.num_x_blocks + i,
                        i * c.cells_x_per_block
                    )
            );
}

template<>
void partition_server::receive_boundary<TOP>(std::size_t step, std::size_t var, future_vector& recv_futures)
{
    if (recv_buffer_top_[var].valid_)
        for (std::size_t i = 0; i < recv_futures.size(); ++i)
            recv_futures[i] =
                hpx::async(
                    hpx::util::bind(
                        boost::ref(recv_buffer_top_[var]),
                        boost::ref(data_[var]),
                        step * c.num_x_blocks + i,
                        i * c.cells_x_per_block
                    )
            );
}

template<>
void partition_server::receive_boundary<BOTTOM_RIGHT>(std::size_t step, std::size_t var, future_vector& recv_futures)
{
    if (recv_buffer_bottom_right_[var].valid_)
        recv_futures[0] =
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
        recv_futures[0] =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_top_left_[var]),
                    boost::ref(data_[var]),
                    step
                )
        );
}

template<std::size_t var>
void partition_server::receive_left_and_bottom_boundaries(std::size_t step, future_grid& recv_futures)
{
    receive_boundary<LEFT>(step, var, recv_futures[LEFT]);
    receive_boundary<BOTTOM>(step, var, recv_futures[BOTTOM]);

    for (uint i = 0; i < recv_futures[RIGHT].size(); ++i)
        recv_futures[RIGHT][i] = hpx::make_ready_future();

    for (uint i = 0; i < recv_futures[TOP].size(); ++i)
        recv_futures[TOP][i] = hpx::make_ready_future();

    recv_futures[BOTTOM_RIGHT][0] = hpx::make_ready_future();
    recv_futures[TOP_LEFT][0] = hpx::make_ready_future();
}

template<std::size_t var>
void partition_server::receive_cross_boundaries(std::size_t step, future_grid& recv_futures)
{
    receive_boundary<LEFT>(step, var, recv_futures[LEFT]);
    receive_boundary<RIGHT>(step, var, recv_futures[RIGHT]);
    receive_boundary<BOTTOM>(step, var, recv_futures[BOTTOM]);
    receive_boundary<TOP>(step, var, recv_futures[TOP]);

    recv_futures[BOTTOM_RIGHT][0] = hpx::make_ready_future();
    recv_futures[TOP_LEFT][0] = hpx::make_ready_future();
}

template<std::size_t var>
void partition_server::receive_all_boundaries(std::size_t step, future_grid& recv_futures)
{
    receive_boundary<LEFT>(step, var, recv_futures[LEFT]);
    receive_boundary<RIGHT>(step, var, recv_futures[RIGHT]);
    receive_boundary<BOTTOM>(step, var, recv_futures[BOTTOM]);
    receive_boundary<TOP>(step, var, recv_futures[TOP]);
    receive_boundary<BOTTOM_RIGHT>(step, var, recv_futures[BOTTOM_RIGHT]);
    receive_boundary<TOP_LEFT>(step, var, recv_futures[TOP_LEFT]);
}

void partition_server::wait_all_boundaries(future_grid& recv_futures)
{
    if (!is_left_)
        hpx::wait_all(recv_futures[LEFT]);

    if (!is_right_)
        hpx::wait_all(recv_futures[RIGHT]);

    if (!is_bottom_)
    {
        hpx::wait_all(recv_futures[BOTTOM]);

        if (!is_right_)
            recv_futures[BOTTOM_RIGHT][0].wait();
    }

    if (!is_top_)
    {
        hpx::wait_all(recv_futures[TOP]);

        if (!is_left_)
            recv_futures[TOP_LEFT][0].wait();
    }
}

template<>
hpx::shared_future<void> partition_server::get_dependency<LEFT>(std::size_t idx_block,
    std::size_t idy_block, future_grid const& recv_futures,
    partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if (is_left_ && idx_block == 0)
        return hpx::make_ready_future();
    else if (idx_block == 0)
        return recv_futures[LEFT][idy_block];

    return calc_futures(idx_block - 1, idy_block);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<RIGHT>(std::size_t idx_block,
    std::size_t idy_block, future_grid const& recv_futures,
    partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if (is_right_ && idx_block == calc_futures.size_x_ - 1)
        return hpx::make_ready_future();
    else if (idx_block == calc_futures.size_x_ - 1)
        return recv_futures[RIGHT][idy_block];

    return calc_futures(idx_block + 1, idy_block);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<BOTTOM>(std::size_t idx_block,
    std::size_t idy_block, future_grid const& recv_futures,
    partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if (is_bottom_ && idy_block == 0)
        return hpx::make_ready_future();

    else if (idy_block == 0)
        return recv_futures[BOTTOM][idx_block];

    return calc_futures(idx_block, idy_block - 1);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<TOP>(std::size_t idx_block,
    std::size_t idy_block, future_grid const& recv_futures,
    partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if (is_top_ && idy_block == calc_futures.size_y_ - 1)
        return hpx::make_ready_future();
    else if (idy_block == calc_futures.size_y_ - 1)
        return recv_futures[TOP][idx_block];

    return calc_futures(idx_block, idy_block + 1);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<BOTTOM_RIGHT>(std::size_t idx_block,
    std::size_t idy_block, future_grid const& recv_futures,
    partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if ((is_bottom_ && idy_block == 0) || (is_right_ && idx_block == calc_futures.size_x_ - 1))
        return hpx::make_ready_future();

    if (idy_block == 0 && idx_block == calc_futures.size_x_ - 1)
        return recv_futures[BOTTOM_RIGHT][0];
    else if (idy_block == 0)
        return recv_futures[BOTTOM][idx_block + 1];
    else if (idx_block == calc_futures.size_x_ - 1)
        return recv_futures[RIGHT][idy_block - 1];

    return calc_futures(idx_block + 1, idy_block - 1);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<TOP_LEFT>(std::size_t idx_block,
    std::size_t idy_block, future_grid const& recv_futures,
    partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if ((is_left_ && idx_block == 0) || (is_top_ && idy_block == calc_futures.size_y_ - 1))
        return hpx::make_ready_future();

    if (idy_block == calc_futures.size_y_ - 1 && idx_block == 0)
        return recv_futures[TOP_LEFT][0];
    else if (idy_block == calc_futures.size_y_ - 1)
        return recv_futures[TOP][idx_block - 1];
    else if (idx_block == 0)
        return recv_futures[LEFT][idy_block + 1];

    return calc_futures(idx_block - 1, idy_block + 1);
}

std::pair<Real, Real> partition_server::do_timestep(Real dt)
{
    hpx::util::high_resolution_timer t;

    clear(send_futures_U);
    clear(send_futures_V);

    for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
    {
        for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
        {
            hpx::shared_future<void> calc_future1 =
                hpx::async(
                    hpx::util::bind(
                        &stencils<STENCIL_SET_VELOCITY_BOUNDARY>::call,
                        boost::ref(data_[U]), boost::ref(data_[V]),
                        boost::ref(data_[U]), boost::ref(data_[V]),
                        boost::ref(cell_type_data_), boost::ref(boundary_cells_(nx_block, ny_block)),
                        c.u_bnd, c.v_bnd, c.bnd_type
                    )
                );

            hpx::shared_future<void> calc_future2 =
                hpx::async(
                    hpx::util::bind(
                        &stencils<STENCIL_SET_VELOCITY_OBSTACLE>::call,
                        boost::ref(data_[U]), boost::ref(data_[V]),
                        boost::ref(data_[U]), boost::ref(data_[V]),
                        boost::ref(cell_type_data_), boost::ref(obstacle_cells_(nx_block, ny_block))
                    )
                );

            hpx::shared_future<void> calc_future3 =
                hpx::async(
                    hpx::util::bind(
                        &stencils<STENCIL_SET_VELOCITY_FLUID>::call,
                        boost::ref(data_[U]), boost::ref(data_[V]),
                        boost::ref(data_[U]), boost::ref(data_[V]),
                        boost::ref(cell_type_data_), boost::ref(fluid_cells_(nx_block, ny_block))
                    )
                );

            hpx::shared_future<void> calc_future =
                hpx::when_all(calc_future1, calc_future2, calc_future3).then(
                    [](hpx::lcos::future<hpx::util::tuple<hpx::lcos::shared_future<void>, hpx::lcos::shared_future<void>, hpx::lcos::shared_future<void> > >)
                    {return;});

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

    receive_all_boundaries<U>(step_, recv_futures_U);
    receive_all_boundaries<V>(step_, recv_futures_V);

    if (c.vtk && next_out_ < t_)
    {
        next_out_ += c.delta_vec;

        hpx::when_all(set_velocity_futures.data_).then(
            hpx::launch::async,
            hpx::util::bind(
                &io::writer::write_vtk,
                boost::ref(data_[P]), boost::ref(data_[U]), boost::ref(data_[V]), boost::ref(cell_type_data_),
                c.num_partitions_x, c.num_partitions_y, c.i_max, c.j_max, c.dx, c.dx, outcount_++,
                c.rank
            )
        ).wait();
    }

    clear(send_futures_F);
    clear(send_futures_G);

    for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
    {
        for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
        {
            hpx::shared_future<void> calc_future1 =
                hpx::dataflow(
                    hpx::util::unwrapped(
                        hpx::util::bind(
                            &stencils<STENCIL_COMPUTE_FG_BOUNDARY_AND_OBSTACLE>::call,
                            boost::ref(data_[F]), boost::ref(data_[G]),
                            boost::ref(data_[U]), boost::ref(data_[V]),
                            boost::ref(cell_type_data_),
                            boost::ref(boundary_cells_(nx_block, ny_block)),
                            boost::ref(obstacle_cells_(nx_block, ny_block))
                        )
                    )
                    , set_velocity_futures(nx_block, ny_block)
                    , get_dependency<LEFT>(nx_block, ny_block, recv_futures_U, set_velocity_futures)
                    , get_dependency<LEFT>(nx_block, ny_block, recv_futures_V, set_velocity_futures)
                    , get_dependency<RIGHT>(nx_block, ny_block, recv_futures_U, set_velocity_futures)
                    , get_dependency<RIGHT>(nx_block, ny_block, recv_futures_V, set_velocity_futures)
                    , get_dependency<BOTTOM>(nx_block, ny_block, recv_futures_U, set_velocity_futures)
                    , get_dependency<BOTTOM>(nx_block, ny_block, recv_futures_V, set_velocity_futures)
                    , get_dependency<TOP>(nx_block, ny_block, recv_futures_U, set_velocity_futures)
                    , get_dependency<TOP>(nx_block, ny_block, recv_futures_V, set_velocity_futures)
                    , get_dependency<TOP_LEFT>(nx_block, ny_block, recv_futures_U, set_velocity_futures)
                    , get_dependency<BOTTOM_RIGHT>(nx_block, ny_block, recv_futures_V, set_velocity_futures)
                );

            hpx::shared_future<void> calc_future2 =
                hpx::dataflow(
                    hpx::util::unwrapped(
                        hpx::util::bind(
                            &stencils<STENCIL_COMPUTE_FG_FLUID>::call,
                            boost::ref(data_[F]), boost::ref(data_[G]),
                            boost::ref(data_[U]), boost::ref(data_[V]),
                            boost::ref(cell_type_data_),
                            boost::ref(fluid_cells_(nx_block, ny_block)),
                            c.re, c.gx, c.gy, c.beta, c.dx, c.dy, dt, c.alpha
                        )
                    )
                    , set_velocity_futures(nx_block, ny_block)
                    , get_dependency<LEFT>(nx_block, ny_block, recv_futures_U, set_velocity_futures)
                    , get_dependency<LEFT>(nx_block, ny_block, recv_futures_V, set_velocity_futures)
                    , get_dependency<RIGHT>(nx_block, ny_block, recv_futures_U, set_velocity_futures)
                    , get_dependency<RIGHT>(nx_block, ny_block, recv_futures_V, set_velocity_futures)
                    , get_dependency<BOTTOM>(nx_block, ny_block, recv_futures_U, set_velocity_futures)
                    , get_dependency<BOTTOM>(nx_block, ny_block, recv_futures_V, set_velocity_futures)
                    , get_dependency<TOP>(nx_block, ny_block, recv_futures_U, set_velocity_futures)
                    , get_dependency<TOP>(nx_block, ny_block, recv_futures_V, set_velocity_futures)
                    , get_dependency<TOP_LEFT>(nx_block, ny_block, recv_futures_U, set_velocity_futures)
                    , get_dependency<BOTTOM_RIGHT>(nx_block, ny_block, recv_futures_V, set_velocity_futures)
                );

            hpx::shared_future<void> calc_future =
                hpx::when_all(calc_future1, calc_future2).then(
                    [](hpx::lcos::future<hpx::util::tuple<hpx::lcos::shared_future<void>, hpx::lcos::shared_future<void>> >)
                    {return;});

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

    receive_left_and_bottom_boundaries<F>(step_, recv_futures_F);
    receive_left_and_bottom_boundaries<G>(step_, recv_futures_G);

    for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
    {
        for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
        {
            compute_rhs_futures(nx_block, ny_block) =
                hpx::dataflow(
                    hpx::util::unwrapped(
                        hpx::util::bind(
                            &stencils<STENCIL_COMPUTE_RHS>::call,
                            boost::ref(rhs_data_), boost::ref(data_[F]),
                            boost::ref(data_[G]), boost::ref(cell_type_data_),
                            boost::ref(fluid_cells_(nx_block, ny_block)),
                            c.dx, c.dy, dt
                        )
                    )
                    , compute_fg_futures(nx_block, ny_block)
                    , get_dependency<LEFT>(nx_block, ny_block, recv_futures_F, compute_fg_futures)
                    , get_dependency<BOTTOM>(nx_block, ny_block, recv_futures_G, compute_fg_futures)
                );
        }
    }

    //TODO make calc_futures members
    hpx::util::high_resolution_timer t1;

    Real rres;
    token.reset();

   // if (!is_left_)
   //     hpx::wait_all(g[LEFT]);

    for (std::size_t iter = 0; iter < c.iter_max; ++iter)
    {
        clear(send_futures_P);

#ifdef WITH_SOR
        for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
        {
            for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
            {
                set_p_futures(nx_block, ny_block) =
                    hpx::dataflow(
                        hpx::util::unwrapped(
                            hpx::util::bind(
                                &stencils<STENCIL_SET_P>::call,
                                boost::ref(data_[P]), boost::ref(cell_type_data_),
                                boost::ref(boundary_cells_(nx_block, ny_block)),
                                boost::ref(obstacle_cells_(nx_block, ny_block)),
                                token
                            )
                        )
                        , sor_cycle_futures[last](nx_block, ny_block)
                        , get_dependency<LEFT>(nx_block, ny_block, recv_futures_P[last], sor_cycle_futures[last])
                        , get_dependency<BOTTOM>(nx_block, ny_block, recv_futures_P[last], sor_cycle_futures[last])
                        , get_dependency<RIGHT>(nx_block, ny_block, recv_futures_P[last], sor_cycle_futures[last])
                        , get_dependency<TOP>(nx_block, ny_block, recv_futures_P[last], sor_cycle_futures[last])
                        , compute_res_futures(nx_block, ny_block).then([](hpx::shared_future<Real> a) {return;})
                    );
            }
        }

        receive_cross_boundaries<P>(iter, recv_futures_P[current]);

        for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
        {
            for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
            {
                hpx::shared_future<void> calc_future =
                    hpx::dataflow(
                        hpx::util::unwrapped(
                            hpx::util::bind(
                                &stencils<STENCIL_SOR>::call,
                                boost::ref(data_[P]), boost::ref(rhs_data_),
                                boost::ref(fluid_cells_(nx_block, ny_block)),
                                c.part1, c.part2, c.dx_sq, c.dy_sq,
                                token, nx_block, ny_block, iter
                            )
                        )
                        , set_p_futures(nx_block, ny_block)
                        , compute_rhs_futures(nx_block, ny_block)
                        , get_dependency<LEFT>(nx_block, ny_block, recv_futures_P[current], sor_cycle_futures[current])
                        , get_dependency<BOTTOM>(nx_block, ny_block, recv_futures_P[current], sor_cycle_futures[current])
                        , get_dependency<RIGHT>(nx_block, ny_block, recv_futures_P[last], sor_cycle_futures[last])
                        , get_dependency<TOP>(nx_block, ny_block, recv_futures_P[last], sor_cycle_futures[last])
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

        send_cross_boundaries<P>(iter, send_futures_P);

        for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
        {
            for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
            {
                compute_res_futures(nx_block, ny_block) =
                    hpx::dataflow(
                        hpx::util::unwrapped(
                            hpx::util::bind(
                                &stencils<STENCIL_COMPUTE_RESIDUAL>::call,
                                boost::ref(data_[P]), boost::ref(rhs_data_),
                                boost::ref(fluid_cells_(nx_block, ny_block)),
                                c.over_dx_sq, c.over_dy_sq,
                                token
                            )
                        )
                        , sor_cycle_futures[current](nx_block, ny_block)
                        , get_dependency<RIGHT>(nx_block, ny_block, recv_futures_P[current], sor_cycle_futures[current])
                        , get_dependency<TOP>(nx_block, ny_block, recv_futures_P[current], sor_cycle_futures[current])
                    );
            }
        }
#else
        for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
        {
            for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
            {
                set_p_futures(nx_block, ny_block) =
                    hpx::dataflow(
                        hpx::util::unwrapped(
                            hpx::util::bind(
                                &stencils<STENCIL_SET_P>::call,
                                boost::ref(data_[P]), boost::ref(cell_type_data_),
                                boost::ref(boundary_cells_(nx_block, ny_block)), boost::ref(obstacle_cells_(nx_block, ny_block)),
                                token
                            )
                        )
                        , sor_cycle_futures[last](nx_block, ny_block)
                        , get_dependency<LEFT>(nx_block, ny_block, recv_futures_P[last], sor_cycle_futures[last])
                        , get_dependency<BOTTOM>(nx_block, ny_block, recv_futures_P[last], sor_cycle_futures[last])
                        , get_dependency<RIGHT>(nx_block, ny_block, recv_futures_P[last], sor_cycle_futures[last])
                        , get_dependency<TOP>(nx_block, ny_block, recv_futures_P[last], sor_cycle_futures[last])
                        , compute_res_futures(nx_block, ny_block).then([](hpx::shared_future<Real> a) {return;})
                    );
            }
        }

        for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
        {
            for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
            {
                hpx::shared_future<void> calc_future =
                    hpx::dataflow(
                        hpx::util::unwrapped(
                            hpx::util::bind(
                                &stencils<STENCIL_JACOBI>::call,
                                boost::ref(data_[P]), boost::ref(rhs_data_),
                                boost::ref(fluid_cells_(nx_block, ny_block)), c.dx_sq, c.dy_sq,
                                token, nx_block, ny_block, iter
                            )
                        )
                        , set_p_futures(nx_block, ny_block)
                        , compute_rhs_futures(nx_block, ny_block)
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

        send_cross_boundaries<P>(iter, send_futures_P);
        receive_cross_boundaries<P>(iter, recv_futures_P[current]);

        for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
        {
            for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
            {
                compute_res_futures(nx_block, ny_block) =
                    hpx::dataflow(
                        hpx::util::unwrapped(
                            hpx::util::bind(
                                &stencils<STENCIL_COMPUTE_RESIDUAL>::call,
                                boost::ref(data_[P]), boost::ref(rhs_data_),
                                boost::ref(fluid_cells_(nx_block, ny_block)), c.over_dx_sq, c.over_dy_sq,
                               token
                            )
                        )
                        , sor_cycle_futures[current](nx_block, ny_block)
                        , get_dependency<LEFT>(nx_block, ny_block, recv_futures_P[current], sor_cycle_futures[current])
                        , get_dependency<RIGHT>(nx_block, ny_block, recv_futures_P[current], sor_cycle_futures[current])
                        , get_dependency<BOTTOM>(nx_block, ny_block, recv_futures_P[current], sor_cycle_futures[current])
                        , get_dependency<TOP>(nx_block, ny_block, recv_futures_P[current], sor_cycle_futures[current])
                    );
            }
        }

#endif
        std::swap(current, last);

        hpx::future<Real> local_residual =
            hpx::dataflow(
                hpx::util::unwrapped(
                    [num_fluid_cells = c.num_fluid_cells](std::vector<Real> residuals)
                    -> Real
                    {
                        Real sum = 0;

                        for (std::size_t i = 0; i < residuals.size(); ++i)
                            sum += residuals[i];

                        return sum / num_fluid_cells;
                    }
                )
                , compute_res_futures.data_
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

                            return std::sqrt(residual);
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
                                        step_ * c.iter_max + iter);
    }

    local_max_uvs.clear();
    local_max_uvs.reserve(c.num_x_blocks * c.num_y_blocks);

    for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
    {
        for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
        {
            local_max_uvs.push_back(
                hpx::dataflow(
                    hpx::util::unwrapped(
                        hpx::util::bind(
                            &stencils<STENCIL_UPDATE_VELOCITY>::call,
                            boost::ref(data_[U]), boost::ref(data_[V]), boost::ref(data_[F]),
                            boost::ref(data_[G]), boost::ref(data_[P]),
                            boost::ref(cell_type_data_),
                            boost::ref(fluid_cells_(nx_block, ny_block)),
                            dt, c.over_dx, c.over_dy
                        )
                    )
                    , get_dependency<LEFT>(nx_block, ny_block, recv_futures_P[last], sor_cycle_futures[last])
                    , get_dependency<RIGHT>(nx_block, ny_block, recv_futures_P[last], sor_cycle_futures[last])
                    , get_dependency<BOTTOM>(nx_block, ny_block, recv_futures_P[last], sor_cycle_futures[last])
                    , get_dependency<TOP>(nx_block, ny_block, recv_futures_P[last], sor_cycle_futures[last])
                    , sor_cycle_futures[last](nx_block, ny_block)
                )
            );
        }
    }

    hpx::future<std::pair<Real, Real> > local_max_uv =
            hpx::dataflow(
                hpx::util::unwrapped(
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
                , local_max_uvs
            );

    t_ += dt;
    ++step_;

    hpx::wait_all(compute_res_futures.data_);

    return local_max_uv.get();
}

}
}
}
