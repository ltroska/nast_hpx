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

HPX_REGISTER_GATHER(double, partition_server_residual_gather);

namespace nast_hpx { namespace grid { namespace server {

partition_server::partition_server(io::config const& cfg)
:   c(cfg),
    cells_x_(c.cells_x_per_partition + 2),
    cells_y_(c.cells_y_per_partition + 2),
    cells_z_(c.cells_z_per_partition + 2),
    is_left_(c.idx == 0),
    is_right_(c.idx == c.num_localities_x - 1),
    is_bottom_(c.idz == 0),
    is_top_(c.idz == c.num_localities_z - 1),
    is_front_(c.idy == 0),
    is_back_(c.idy == c.num_localities_y - 1)
{
    if (c.verbose && c.rank == 0)
    {
        std::cout << "Solver: blockwise Jacobi" << std::endl;
        std::cout << "Parellelization: custom grain size" << std::endl;
        std::cout << c << std::endl;
    }

    step_ = 0;

    data_[U].resize(cells_x_, cells_y_, cells_z_, 0);
    data_[V].resize(cells_x_, cells_y_, cells_z_, 0);
    data_[W].resize(cells_x_, cells_y_, cells_z_, 0);
    data_[F].resize(cells_x_, cells_y_, cells_z_, 0);
    data_[G].resize(cells_x_, cells_y_, cells_z_, 0);
    data_[H].resize(cells_x_, cells_y_, cells_z_, 0);
    data_[P].resize(cells_x_, cells_y_, cells_z_, 0);
    rhs_data_.resize(cells_x_, cells_y_, cells_z_, 0);

    cell_type_data_.resize(cells_x_, cells_y_, cells_z_);

    for (std::size_t k = 1; k < cells_z_ - 1; ++k)
        for (std::size_t j = 1; j < cells_y_ - 1; ++j)
            for (std::size_t i = 1; i < cells_x_ - 1; ++i)
            {
                cell_type_data_(i, j, k) = std::move(c.flag_grid[k * cells_x_ * cells_y_ + j * cells_x_ + i]);

                if (cell_type_data_(i, j, k)[is_fluid])
                    fluid_cells_.emplace_back(i, j, k);
                else if ( (cell_type_data_(i, j, k)[is_obstacle] && !cell_type_data_(i, j, k)[is_boundary] && cell_type_data_(i,  j, k).count() > 1)
                            || cell_type_data_(i, j, k).count() > 2)
                    obstacle_cells_.emplace_back(i, j, k);
            }

    fluid_stride = fluid_cells_.size() / c.threads + (fluid_cells_.size() % c.threads > 0);
    obstacle_stride = obstacle_cells_.size() / c.threads + (obstacle_cells_.size() % c.threads > 0);
}

std::size_t partition_server::get_id(std::size_t dir_x, std::size_t dir_y, std::size_t dir_z)
{
    return (c.idz + dir_z) * c.num_localities_x * c.num_localities_y + (c.idy + dir_y) * c.num_localities_x + (c.idx + dir_x);
}

void partition_server::init()
{
    for (std::size_t var = 0; var < NUM_VARIABLES; ++var)
        data_[var].clear(0);

    rhs_data_.clear(0);

    step_ = 0;
    t_ = 0;
    next_out_ = -1e-12;
    outcount_ = 0;

    std::vector<hpx::future<hpx::id_type > > parts =
        hpx::find_all_from_basename(partition_basename, c.num_localities);

    ids_ = hpx::when_all(parts).then(hpx::util::unwrapped2(
                [](std::vector<hpx::id_type>&& ids) -> std::vector<hpx::id_type>
                { return ids;})
            ).get();

    if (!is_left_)
    {
        send_buffer_left_.dest_ = hpx::find_from_basename(partition_basename, get_id(-1, 0, 0)).get();
        recv_buffer_left_[U].valid_ = true;
        recv_buffer_left_[V].valid_ = true;
        recv_buffer_left_[W].valid_ = true;
        recv_buffer_left_[F].valid_ = true;
        recv_buffer_left_[P].valid_ = true;
    }

    if (!is_right_)
    {
        send_buffer_right_.dest_ = hpx::find_from_basename(partition_basename, get_id(1, 0, 0)).get();

        recv_buffer_right_[U].valid_ = true;
        recv_buffer_right_[V].valid_ = true;
        recv_buffer_right_[W].valid_ = true;
        recv_buffer_right_[P].valid_ = true;
    }

    if (!is_bottom_)
    {
        send_buffer_bottom_.dest_ = hpx::find_from_basename(partition_basename, get_id(0, 0, -1)).get();

        recv_buffer_bottom_[U].valid_ = true;
        recv_buffer_bottom_[V].valid_ = true;
        recv_buffer_bottom_[W].valid_ = true;
        recv_buffer_bottom_[H].valid_ = true;
        recv_buffer_bottom_[P].valid_ = true;    }

    if (!is_top_)
    {
        send_buffer_top_.dest_ = hpx::find_from_basename(partition_basename, get_id(0, 0, 1)).get();

        recv_buffer_top_[U].valid_ = true;
        recv_buffer_top_[V].valid_ = true;
        recv_buffer_top_[W].valid_ = true;
        recv_buffer_top_[P].valid_ = true;    }

    if (!is_front_)
    {
        send_buffer_front_.dest_ = hpx::find_from_basename(partition_basename, get_id(0, -1, 0)).get();

        recv_buffer_front_[U].valid_ = true;
        recv_buffer_front_[V].valid_ = true;
        recv_buffer_front_[W].valid_ = true;
        recv_buffer_front_[G].valid_ = true;
        recv_buffer_front_[P].valid_ = true;    }

    if (!is_back_)
    {
        send_buffer_back_.dest_ = hpx::find_from_basename(partition_basename, get_id(0, 1, 0)).get();

        recv_buffer_back_[U].valid_ = true;
        recv_buffer_back_[V].valid_ = true;
        recv_buffer_back_[W].valid_ = true;
        recv_buffer_back_[P].valid_ = true;    }

    if (!is_back_ && !is_left_)
    {
        send_buffer_back_left_.dest_ = hpx::find_from_basename(partition_basename, get_id(-1, 1, 0)).get();

        recv_buffer_back_left_[U].valid_ = true;
    }

    if (!is_front_ && !is_right_)
    {
        send_buffer_front_right_.dest_ = hpx::find_from_basename(partition_basename, get_id(1, -1, 0)).get();;

        recv_buffer_front_right_[V].valid_ = true;
    }

    if (!is_bottom_ && !is_right_)
    {
        send_buffer_bottom_right_.dest_ = hpx::find_from_basename(partition_basename, get_id(1, 0, -1)).get();

        recv_buffer_bottom_right_[W].valid_ = true;
    }

    if (!is_top_ && !is_left_)
    {
        send_buffer_top_left_.dest_ = hpx::find_from_basename(partition_basename, get_id(-1, 0, 1)).get();

        recv_buffer_top_left_[U].valid_ = true;
    }

    if (!is_back_ && !is_bottom_)
    {
        send_buffer_back_bottom_.dest_ = hpx::find_from_basename(partition_basename, get_id(0, 1, -1)).get();

        recv_buffer_back_bottom_[W].valid_ = true;
    }

    if (!is_front_ && !is_top_)
    {
        send_buffer_front_top_.dest_ = hpx::find_from_basename(partition_basename, get_id(0, -1, 1)).get();

        recv_buffer_front_top_[V].valid_ = true;
    }

    recv_futures.resize(NUM_VARIABLES);

    for (std::size_t var = 0; var < NUM_VARIABLES; ++var)
        recv_futures[var].resize(NUM_DIRECTIONS);

    set_velocity_futures.resize(c.threads);
    compute_fg_futures.resize(c.threads);
    compute_rhs_futures.resize(c.threads);
    set_p_futures.resize(c.threads);

    compute_res_futures.resize(c.threads);
    for (auto& a : compute_res_futures)
        a = hpx::make_ready_future(0.);

    solver_cycle_futures.resize(c.threads);
    for (auto& a : solver_cycle_futures)
        a = hpx::make_ready_future();

    local_max_uvs.resize(c.threads);

    token.reset();
}

template<>
void partition_server::send_boundary<LEFT>(std::size_t step, std::size_t var, future_vector& send_future)
{
        hpx::when_all(send_future).then(
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
void partition_server::send_boundary<RIGHT>(std::size_t step, std::size_t var, future_vector& send_future)
{
        hpx::when_all(send_future).then(
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
void partition_server::send_boundary<BOTTOM>(std::size_t step, std::size_t var, future_vector& send_future)
{
        hpx::when_all(send_future).then(
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
void partition_server::send_boundary<TOP>(std::size_t step, std::size_t var, future_vector& send_future)
{
        hpx::when_all(send_future).then(
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
void partition_server::send_boundary<FRONT>(std::size_t step, std::size_t var, future_vector& send_future)
{
        hpx::when_all(send_future).then(
            hpx::launch::async,
            hpx::util::bind(
                boost::ref(send_buffer_front_),
                boost::ref(data_[var]),
                step,
                var
            )
        );
}

template<>
void partition_server::send_boundary<BACK>(std::size_t step, std::size_t var, future_vector& send_future)
{
        hpx::when_all(send_future).then(
            hpx::launch::async,
            hpx::util::bind(
                boost::ref(send_buffer_back_),
                boost::ref(data_[var]),
                step,
                var
            )
        );
}

template<>
void partition_server::send_boundary<BACK_LEFT>(std::size_t step, std::size_t var, future_vector& send_future)
{
    hpx::when_all(send_future).then(
        hpx::launch::async,
        hpx::util::bind(
            boost::ref(send_buffer_back_left_),
            boost::ref(data_[var]),
            step,
            var
        )
    );
}

template<>
void partition_server::send_boundary<FRONT_RIGHT>(std::size_t step, std::size_t var, future_vector& send_future)
{
    hpx::when_all(send_future).then(
        hpx::launch::async,
        hpx::util::bind(
            boost::ref(send_buffer_front_right_),
            boost::ref(data_[var]),
            step,
            var
        )
    );
}

template<>
void partition_server::send_boundary<BOTTOM_RIGHT>(std::size_t step, std::size_t var, future_vector& send_future)
{
    hpx::when_all(send_future).then(
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
void partition_server::send_boundary<TOP_LEFT>(std::size_t step, std::size_t var, future_vector& send_future)
{
    hpx::when_all(send_future).then(
        hpx::launch::async,
        hpx::util::bind(
            boost::ref(send_buffer_top_left_),
            boost::ref(data_[var]),
            step,
            var
        )
    );
}

template<>
void partition_server::send_boundary<BACK_BOTTOM>(std::size_t step, std::size_t var,future_vector& send_future)
{
    hpx::when_all(send_future).then(
        hpx::launch::async,
        hpx::util::bind(
            boost::ref(send_buffer_back_bottom_),
            boost::ref(data_[var]),
            step,
            var
        )
    );
}

template<>
void partition_server::send_boundary<FRONT_TOP>(std::size_t step, std::size_t var, future_vector& send_future)
{
    hpx::when_all(send_future).then(
        hpx::launch::async,
        hpx::util::bind(
            boost::ref(send_buffer_front_top_),
            boost::ref(data_[var]),
            step,
            var
        )
    );
}

template<>
void partition_server::receive_boundary<LEFT>(std::size_t step, std::size_t var, future_grid& recv_futures)
{
        recv_futures[var][LEFT] =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_left_[var]),
                    boost::ref(data_[var]),
                    step
                )
            );
}

template<>
void partition_server::receive_boundary<RIGHT>(std::size_t step, std::size_t var, future_grid& recv_futures)
{
        recv_futures[var][RIGHT] =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_right_[var]),
                    boost::ref(data_[var]),
                    step
                )
            );
}

template<>
void partition_server::receive_boundary<BOTTOM>(std::size_t step, std::size_t var, future_grid& recv_futures)
{
        recv_futures[var][BOTTOM] =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_bottom_[var]),
                    boost::ref(data_[var]),
                    step
                )
            );
}

template<>
void partition_server::receive_boundary<TOP>(std::size_t step, std::size_t var, future_grid& recv_futures)
{
        recv_futures[var][TOP] =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_top_[var]),
                    boost::ref(data_[var]),
                    step
                )
            );
}

template<>
void partition_server::receive_boundary<FRONT>(std::size_t step, std::size_t var, future_grid& recv_futures)
{
        recv_futures[var][FRONT] =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_front_[var]),
                    boost::ref(data_[var]),
                    step
                )
            );
}

template<>
void partition_server::receive_boundary<BACK>(std::size_t step, std::size_t var, future_grid& recv_futures)
{
        recv_futures[var][BACK] =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_back_[var]),
                    boost::ref(data_[var]),
                    step
                )
            );
}

template<>
void partition_server::receive_boundary<BACK_LEFT>(std::size_t step, std::size_t var, future_grid& recv_futures)
{
        recv_futures[var][BACK_LEFT] =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_back_left_[var]),
                    boost::ref(data_[var]),
                    step
                )
            );
}

template<>
void partition_server::receive_boundary<FRONT_RIGHT>(std::size_t step, std::size_t var, future_grid& recv_futures)
{
        recv_futures[var][FRONT_RIGHT] =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_front_right_[var]),
                    boost::ref(data_[var]),
                    step
                )
            );
}

template<>
void partition_server::receive_boundary<BOTTOM_RIGHT>(std::size_t step, std::size_t var, future_grid& recv_futures)
{
        recv_futures[var][BOTTOM_RIGHT] =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_bottom_right_[var]),
                    boost::ref(data_[var]),
                    step
                )
            );
}

template<>
void partition_server::receive_boundary<TOP_LEFT>(std::size_t step, std::size_t var, future_grid& recv_futures)
{
    recv_futures[var][TOP_LEFT] =
        hpx::async(
            hpx::util::bind(
                boost::ref(recv_buffer_top_left_[var]),
                boost::ref(data_[var]),
                step
            )
        );
}

template<>
void partition_server::receive_boundary<BACK_BOTTOM>(std::size_t step, std::size_t var, future_grid& recv_futures)
{
    recv_futures[var][BACK_BOTTOM] =
        hpx::async(
            hpx::util::bind(
                boost::ref(recv_buffer_back_bottom_[var]),
                boost::ref(data_[var]),
                step
            )
        );
}

template<>
void partition_server::receive_boundary<FRONT_TOP>(std::size_t step, std::size_t var, future_grid& recv_futures)
{
        recv_futures[var][FRONT_TOP] =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_front_top_[var]),
                    boost::ref(data_[var]),
                    step
                )
            );
}

template<>
hpx::shared_future<void> partition_server::get_dependency<LEFT>(future_vector const& recv_futures)
{
    if (is_left_)
        return hpx::make_ready_future();
    else
        return recv_futures[LEFT];
}

template<>
hpx::shared_future<void> partition_server::get_dependency<RIGHT>(future_vector const& recv_futures)
{
    if (is_right_)
        return hpx::make_ready_future();
    else
        return recv_futures[RIGHT];
}

template<>
hpx::shared_future<void> partition_server::get_dependency<BOTTOM>(future_vector const& recv_futures)
{
    if (is_bottom_)
        return hpx::make_ready_future();
    else
        return recv_futures[BOTTOM];
}

template<>
hpx::shared_future<void> partition_server::get_dependency<TOP>(future_vector const& recv_futures)
{
    if (is_top_)
        return hpx::make_ready_future();
    else
        return recv_futures[TOP];
}

template<>
hpx::shared_future<void> partition_server::get_dependency<FRONT>(future_vector const& recv_futures)
{
    if (is_front_)
        return hpx::make_ready_future();
    else
        return recv_futures[FRONT];
}

template<>
hpx::shared_future<void> partition_server::get_dependency<BACK>(future_vector const& recv_futures)
{
    if (is_back_)
        return hpx::make_ready_future();
    else
        return recv_futures[BACK];
}

template<>
hpx::shared_future<void> partition_server::get_dependency<BACK_LEFT>(future_vector const& recv_futures)
{
    if (is_back_ || is_left_)
        return hpx::make_ready_future();
    else
        return recv_futures[BACK_LEFT];
}

template<>
hpx::shared_future<void> partition_server::get_dependency<FRONT_RIGHT>(future_vector const& recv_futures)
{
    if (is_front_ || is_right_)
        return hpx::make_ready_future();
    else
        return recv_futures[FRONT_RIGHT];
}

template<>
hpx::shared_future<void> partition_server::get_dependency<BOTTOM_RIGHT>(future_vector const& recv_futures)
{
    if (is_bottom_ || is_right_)
        return hpx::make_ready_future();
    else
        return recv_futures[BOTTOM_RIGHT];
}

template<>
hpx::shared_future<void> partition_server::get_dependency<TOP_LEFT>(future_vector const& recv_futures)
{
    if (is_top_ || is_left_)
        return hpx::make_ready_future();
    else
        return recv_futures[TOP_LEFT];
}

template<>
hpx::shared_future<void> partition_server::get_dependency<BACK_BOTTOM>(future_vector const& recv_futures)
{
    if (is_back_ || is_bottom_)
        return hpx::make_ready_future();
    else
        return recv_futures[BACK_BOTTOM];
}

template<>
hpx::shared_future<void> partition_server::get_dependency<FRONT_TOP>(future_vector const& recv_futures)
{
    if (is_front_ || is_top_)
        return hpx::make_ready_future();
    else
        return recv_futures[FRONT_TOP];
}

void partition_server::send_boundaries_U(future_vector& send_future, std::size_t step)
{
    if (!is_left_)
        send_boundary<LEFT>(step, U, send_future);

    if (!is_right_)
        send_boundary<RIGHT>(step, U, send_future);

    if (!is_bottom_)
        send_boundary<BOTTOM>(step, U, send_future);

    if (!is_top_)
        send_boundary<TOP>(step, U, send_future);

    if (!is_front_)
        send_boundary<FRONT>(step, U, send_future);

    if (!is_back_)
        send_boundary<BACK>(step, U, send_future);

    if (!is_front_ && !is_right_)
        send_boundary<FRONT_RIGHT>(step, U, send_future);

    if (!is_bottom_ && !is_right_)
        send_boundary<BOTTOM_RIGHT>(step, U, send_future);
}

void partition_server::receive_boundaries_U(future_grid& recv_futures, std::size_t step)
{
    if (!is_left_)
        receive_boundary<LEFT>(step, U, recv_futures);

    if (!is_right_)
        receive_boundary<RIGHT>(step, U, recv_futures);

    if (!is_bottom_)
        receive_boundary<BOTTOM>(step, U, recv_futures);

    if (!is_top_)
        receive_boundary<TOP>(step, U, recv_futures);

    if (!is_front_)
        receive_boundary<FRONT>(step, U, recv_futures);

    if (!is_back_)
        receive_boundary<BACK>(step, U, recv_futures);

    if (!is_back_ && !is_left_)
        receive_boundary<BACK_LEFT>(step, U, recv_futures);

    if (!is_top_ && !is_left_)
        receive_boundary<TOP_LEFT>(step, U, recv_futures);

}

void partition_server::send_boundaries_V(future_vector& send_future, std::size_t step)
{
    if (!is_left_)
        send_boundary<LEFT>(step, V, send_future);

    if (!is_right_)
        send_boundary<RIGHT>(step, V, send_future);

    if (!is_bottom_)
        send_boundary<BOTTOM>(step, V, send_future);

    if (!is_top_)
        send_boundary<TOP>(step, V, send_future);

    if (!is_front_)
        send_boundary<FRONT>(step, V, send_future);

    if (!is_back_)
        send_boundary<BACK>(step, V, send_future);

    if (!is_back_ && !is_left_)
        send_boundary<BACK_LEFT>(step, V, send_future);

    if (!is_back_ && !is_bottom_)
        send_boundary<BACK_BOTTOM>(step, V, send_future);
}

void partition_server::receive_boundaries_V(future_grid& recv_futures, std::size_t step)
{
    if (!is_left_)
        receive_boundary<LEFT>(step, V, recv_futures);

    if (!is_right_)
        receive_boundary<RIGHT>(step, V, recv_futures);

    if (!is_bottom_)
        receive_boundary<BOTTOM>(step, V, recv_futures);

    if (!is_top_)
        receive_boundary<TOP>(step, V, recv_futures);

    if (!is_front_)
        receive_boundary<FRONT>(step, V, recv_futures);

    if (!is_back_)
        receive_boundary<BACK>(step, V, recv_futures);

    if (!is_front_ && !is_right_)
        receive_boundary<FRONT_RIGHT>(step, V, recv_futures);

    if (!is_front_ && !is_top_)
        receive_boundary<FRONT_TOP>(step, V, recv_futures);
}

void partition_server::send_boundaries_W(future_vector& send_future, std::size_t step)
{
    if (!is_left_)
        send_boundary<LEFT>(step, W, send_future);

    if (!is_right_)
        send_boundary<RIGHT>(step, W, send_future);

    if (!is_bottom_)
        send_boundary<BOTTOM>(step, W, send_future);

    if (!is_top_)
        send_boundary<TOP>(step, W, send_future);

    if (!is_front_)
        send_boundary<FRONT>(step, W, send_future);

    if (!is_back_)
        send_boundary<BACK>(step, W, send_future);

    if (!is_top_ && !is_left_)
        send_boundary<TOP_LEFT>(step, W, send_future);

    if (!is_top_ && !is_front_)
        send_boundary<FRONT_TOP>(step, W, send_future);
}

void partition_server::receive_boundaries_W(future_grid& recv_futures, std::size_t step)
{
    if (!is_left_)
        receive_boundary<LEFT>(step, W, recv_futures);

    if (!is_right_)
        receive_boundary<RIGHT>(step, W, recv_futures);

    if (!is_bottom_)
        receive_boundary<BOTTOM>(step, W, recv_futures);

    if (!is_top_)
        receive_boundary<TOP>(step, W, recv_futures);

    if (!is_front_)
        receive_boundary<FRONT>(step, W, recv_futures);

    if (!is_back_)
        receive_boundary<BACK>(step, W, recv_futures);

    if (!is_bottom_ && !is_right_)
        receive_boundary<BOTTOM_RIGHT>(step, W, recv_futures);

    if (!is_bottom_ && !is_back_)
        receive_boundary<BACK_BOTTOM>(step, W, recv_futures);
}

void partition_server::send_boundaries_F(future_vector& send_future, std::size_t step)
{
    if (!is_right_)
        send_boundary<RIGHT>(step, F, send_future);
}

void partition_server::receive_boundaries_F(future_grid& recv_futures, std::size_t step)
{
    if (!is_left_)
        receive_boundary<LEFT>(step, F, recv_futures);
}

void partition_server::send_boundaries_G(future_vector& send_future, std::size_t step)
{
    if (!is_back_)
        send_boundary<BACK>(step, G, send_future);
}

void partition_server::receive_boundaries_G(future_grid& recv_futures, std::size_t step)
{
    if (!is_front_)
        receive_boundary<FRONT>(step, G, recv_futures);
}

void partition_server::send_boundaries_H(future_vector& send_future, std::size_t step)
{
    if (!is_top_)
        send_boundary<TOP>(step, H, send_future);
}


void partition_server::receive_boundaries_H(future_grid& recv_futures, std::size_t step)
{
    if (!is_bottom_)
        receive_boundary<BOTTOM>(step, H, recv_futures);
}

void partition_server::send_boundaries_P(future_vector& send_future, std::size_t step)
{
    if (!is_left_)
        send_boundary<LEFT>(step, P, send_future);

    if (!is_right_)
        send_boundary<RIGHT>(step, P, send_future);

    if (!is_bottom_)
        send_boundary<BOTTOM>(step, P, send_future);

    if (!is_top_)
        send_boundary<TOP>(step, P, send_future);

    if (!is_front_)
        send_boundary<FRONT>(step, P, send_future);

    if (!is_back_)
        send_boundary<BACK>(step, P, send_future);
}

void partition_server::receive_boundaries_P(future_grid& recv_futures, std::size_t step)
{
    if (!is_left_)
        receive_boundary<LEFT>(step, P, recv_futures);

    if (!is_right_)
        receive_boundary<RIGHT>(step, P, recv_futures);

    if (!is_bottom_)
        receive_boundary<BOTTOM>(step, P, recv_futures);

    if (!is_top_)
        receive_boundary<TOP>(step, P, recv_futures);

    if (!is_front_)
        receive_boundary<FRONT>(step, P, recv_futures);

    if (!is_back_)
        receive_boundary<BACK>(step, P, recv_futures);
}

template<typename Iter> inline
Iter safe_advance(Iter it, Iter end, std::size_t stride)
{
    return (stride > static_cast<std::size_t>(end - it)) ? end : it + stride;
}

hpx::future<triple<double> > partition_server::do_timestep(double dt)
{
    auto beginObstacle = obstacle_cells_.begin();
    auto endObstacle = safe_advance(beginObstacle, obstacle_cells_.end(), obstacle_stride);

    for (std::size_t thread = 0; thread < c.threads; ++thread)
    {
        set_velocity_futures[thread] =
            hpx::async(
                hpx::util::bind(
                    &stencils<STENCIL_SET_VELOCITY_OBSTACLE>::call,
                    boost::ref(data_[U]), boost::ref(data_[V]), boost::ref(data_[W]),
                    boost::ref(cell_type_data_),
                    beginObstacle, endObstacle,
                    boost::ref(c.bnd_condition)
                )
            );

        beginObstacle = safe_advance(beginObstacle, obstacle_cells_.end(), obstacle_stride);
        endObstacle = safe_advance(endObstacle, obstacle_cells_.end(), obstacle_stride);
    }

    send_boundaries_U(set_velocity_futures, step_);
    receive_boundaries_U(recv_futures, step_);

    send_boundaries_V(set_velocity_futures, step_);
    receive_boundaries_V(recv_futures, step_);

    send_boundaries_W(set_velocity_futures, step_);
    receive_boundaries_W(recv_futures, step_);

    if (c.vtk && next_out_ < t_)
    {
        next_out_ += c.delta_vec;

        if (c.verbose && c.rank == 0)
            std::cout << "Output to .vtk in step " << step_ << std::endl;

        hpx::when_all(set_velocity_futures).then(
            hpx::launch::async,
            hpx::util::bind(
                &io::writer::write_vtk,
                boost::ref(data_[P]), boost::ref(data_[U]), boost::ref(data_[V]), boost::ref(data_[W]), boost::ref(cell_type_data_),
                c.num_localities_x, c.num_localities_y, c.num_localities_z, c.i_max, c.j_max, c.k_max, c.dx, c.dx, c.dz, outcount_++,
                c.rank, c.idx, c.idy, c.idz
            )
        ).wait();
    }

    auto beginFluid = fluid_cells_.begin();
    auto endFluid = safe_advance(beginFluid, fluid_cells_.end(), fluid_stride);
    beginObstacle = obstacle_cells_.begin();
    endObstacle = safe_advance(beginObstacle, obstacle_cells_.end(), obstacle_stride);

    for (std::size_t thread = 0; thread < c.threads; ++thread)
    {
        compute_fg_futures[thread] =
            hpx::dataflow(
                hpx::util::unwrapped(
                    hpx::util::bind(
                        &stencils<STENCIL_COMPUTE_FG>::call,
                        boost::ref(data_[F]), boost::ref(data_[G]), boost::ref(data_[H]),
                        boost::ref(data_[U]), boost::ref(data_[V]), boost::ref(data_[W]),
                        boost::ref(cell_type_data_),
                        beginObstacle, endObstacle,
                        beginFluid, endFluid,
                        c.re, c.gx, c.gy, c.gz, c.dx, c.dy, c.dz,
                        c.dx_sq, c.dy_sq, c.dz_sq, dt, c.alpha
                    )
                )
                , static_cast<hpx::future<void> >(hpx::when_all(set_velocity_futures))
                , get_dependency<LEFT>(recv_futures[U])
                , get_dependency<RIGHT>(recv_futures[U])
                , get_dependency<BOTTOM>(recv_futures[U])
                , get_dependency<TOP>(recv_futures[U])
                , get_dependency<FRONT>(recv_futures[U])
                , get_dependency<BACK>(recv_futures[U])
                , get_dependency<BACK_LEFT>(recv_futures[U])
                , get_dependency<TOP_LEFT>(recv_futures[U])
                , get_dependency<LEFT>(recv_futures[V])
                , get_dependency<RIGHT>(recv_futures[V])
                , get_dependency<BOTTOM>(recv_futures[V])
                , get_dependency<TOP>(recv_futures[V])
                , get_dependency<FRONT>(recv_futures[V])
                , get_dependency<BACK>(recv_futures[V])
                , get_dependency<FRONT_RIGHT>(recv_futures[V])
                , get_dependency<FRONT_TOP>(recv_futures[V])
                , get_dependency<LEFT>(recv_futures[W])
                , get_dependency<RIGHT>(recv_futures[W])
                , get_dependency<BOTTOM>(recv_futures[W])
                , get_dependency<TOP>(recv_futures[W])
                , get_dependency<FRONT>(recv_futures[W])
                , get_dependency<BACK>(recv_futures[W])
                , get_dependency<BOTTOM_RIGHT>(recv_futures[W])
                , get_dependency<BACK_BOTTOM>(recv_futures[W])
            );

        beginFluid = safe_advance(beginFluid, fluid_cells_.end(), fluid_stride);
        endFluid = safe_advance(endFluid, fluid_cells_.end(), fluid_stride);
        beginObstacle = safe_advance(beginObstacle, obstacle_cells_.end(), obstacle_stride);
        endObstacle = safe_advance(endObstacle, obstacle_cells_.end(), obstacle_stride);
    }

    send_boundaries_F(compute_fg_futures, step_);
    receive_boundaries_F(recv_futures, step_);

    send_boundaries_G(compute_fg_futures, step_);
    receive_boundaries_G(recv_futures, step_);

    send_boundaries_H(compute_fg_futures, step_);
    receive_boundaries_H(recv_futures, step_);


    beginFluid = fluid_cells_.begin();
    endFluid = safe_advance(beginFluid, fluid_cells_.end(), fluid_stride);

    for (std::size_t thread = 0; thread < c.threads; ++thread)
    {
        compute_rhs_futures[thread] =
                hpx::dataflow(
                    hpx::util::unwrapped(
                        hpx::util::bind(
                            &stencils<STENCIL_COMPUTE_RHS>::call,
                            boost::ref(rhs_data_),
                            boost::ref(data_[F]), boost::ref(data_[G]), boost::ref(data_[H]),
                            beginFluid, endFluid,
                            c.dx, c.dy, c.dz, dt
                        )
                    )
                    , static_cast<hpx::future<void> >(hpx::when_all(compute_fg_futures))
                    , get_dependency<LEFT>(recv_futures[F])
                    , get_dependency<FRONT>(recv_futures[G])
                    , get_dependency<BOTTOM>(recv_futures[H])
                );

        beginFluid = safe_advance(beginFluid, fluid_cells_.end(), fluid_stride);
        endFluid = safe_advance(endFluid, fluid_cells_.end(), fluid_stride);
    }

    token.reset();
    for (std::size_t iter = 0; iter < c.iter_max; ++iter)
    {
        beginObstacle = obstacle_cells_.begin();
        endObstacle = safe_advance(beginObstacle, obstacle_cells_.end(), obstacle_stride);

        for (std::size_t thread = 0; thread < c.threads; ++thread)
        {
            set_p_futures[thread] =
                hpx::dataflow(
                    hpx::util::unwrapped(
                        hpx::util::bind(
                            &stencils<STENCIL_SET_P_OBSTACLE>::call,
                            boost::ref(data_[P]),
                            boost::ref(cell_type_data_),
                            beginObstacle, endObstacle,
                            token
                        )
                    )
                    , static_cast<hpx::future<void> >(hpx::when_all(compute_rhs_futures))
                    , static_cast<hpx::future<void> >(hpx::when_all(compute_res_futures))
                );

            beginObstacle = safe_advance(beginObstacle, obstacle_cells_.end(), obstacle_stride);
            endObstacle = safe_advance(endObstacle, obstacle_cells_.end(), obstacle_stride);
        }


        beginFluid = fluid_cells_.begin();
        endFluid = safe_advance(beginFluid, fluid_cells_.end(), fluid_stride);

        for (std::size_t thread = 0; thread < c.threads; ++thread)
        {
            solver_cycle_futures[thread] =
                hpx::dataflow(
                    hpx::util::unwrapped(
                        hpx::util::bind(
                            &stencils<STENCIL_JACOBI>::call,
                            boost::ref(data_[P]),
                            boost::ref(rhs_data_),
                            beginFluid, endFluid,
                            c.dx_sq, c.dy_sq, c.dz_sq, token
                        )
                    )
                    , static_cast<hpx::future<void> >(hpx::when_all(set_p_futures))
                );

            beginFluid = safe_advance(beginFluid, fluid_cells_.end(), fluid_stride);
            endFluid = safe_advance(endFluid, fluid_cells_.end(), fluid_stride);
        }

        send_boundaries_P(solver_cycle_futures, step_ * c.iter_max + iter);
        receive_boundaries_P(recv_futures, step_ * c.iter_max + iter);

        beginFluid = fluid_cells_.begin();
        endFluid = safe_advance(beginFluid, fluid_cells_.end(), fluid_stride);

        for (std::size_t thread = 0; thread < c.threads; ++thread)
        {
            compute_res_futures[thread] =
                hpx::dataflow(
                    hpx::util::unwrapped(
                            hpx::util::bind(
                                &stencils<STENCIL_COMPUTE_RESIDUAL>::call,
                                boost::ref(data_[P]),
                                boost::ref(rhs_data_),
                                beginFluid, endFluid,
                                c.dx_sq, c.dy_sq, c.dz_sq, token
                            )
                        )
                        , static_cast<hpx::future<void> >(hpx::when_all(solver_cycle_futures))
                        , get_dependency<LEFT>(recv_futures[P])
                        , get_dependency<RIGHT>(recv_futures[P])
                        , get_dependency<BOTTOM>(recv_futures[P])
                        , get_dependency<TOP>(recv_futures[P])
                        , get_dependency<FRONT>(recv_futures[P])
                        , get_dependency<BACK>(recv_futures[P])
                    );

            beginFluid = safe_advance(beginFluid, fluid_cells_.end(), fluid_stride);
            endFluid = safe_advance(endFluid, fluid_cells_.end(), fluid_stride);
        }

        hpx::future<double> local_residual =
            hpx::dataflow(
                hpx::util::unwrapped(
                    [num_fluid_cells = c.num_fluid_cells](std::vector<double> residuals)
                    -> double
                    {
                        double sum = 0;

                        for (std::size_t i = 0; i < residuals.size(); ++i)
                            sum += residuals[i];

                        return sum / num_fluid_cells;
                    }
                )
                , compute_res_futures
            );

            if (c.rank == 0)
            {
                hpx::future<std::vector<double> > partial_residuals =
                    hpx::lcos::gather_here(residual_basename,
                                            std::move(local_residual),
                                            c.num_localities, step_ * c.iter_max + iter, 0);

                partial_residuals.then(
                    hpx::util::unwrapped(
                        [dt, iter_int = iter, step = step_, t = t_, this](std::vector<double> local_residuals)
                        {
                            double residual = 0;

                            for (std::size_t i = 0; i < local_residuals.size(); ++i)
                                residual += local_residuals[i];

                            residual = std::sqrt(residual);

                            if ((residual < c.eps || iter_int == c.iter_max - 1)
                                && !token.was_cancelled())
                            {
                                if (c.verbose)
                                    std::cout << "step = " << step
                                        << ", t = " << t
                                        << ", dt = " << dt
                                        << ", iter = "<< iter_int
                                        << ", residual = " << residual
                                        << std::endl;

                                hpx::lcos::broadcast_apply<cancel_action>(ids_);
                                token.cancel();
                            }

                        }
                    )
                );
            }
            // if not root locality, send residual to root locality
            else
                hpx::lcos::gather_there(residual_basename, std::move(local_residual),
                                            step_ * c.iter_max + iter, 0, c.rank);
    }

    beginFluid = fluid_cells_.begin();
    endFluid = safe_advance(beginFluid, fluid_cells_.end(), fluid_stride);

    for (std::size_t thread = 0; thread < c.threads; ++thread)
    {
       local_max_uvs[thread] =
            hpx::dataflow(
                hpx::util::unwrapped(
                    hpx::util::bind(
                        &stencils<STENCIL_UPDATE_VELOCITY>::call,
                        boost::ref(data_[U]), boost::ref(data_[V]), boost::ref(data_[W]),
                        boost::ref(data_[F]), boost::ref(data_[G]), boost::ref(data_[H]),
                        boost::ref(data_[P]),
                        boost::ref(cell_type_data_),
                        beginFluid, endFluid,
                        dt, c.over_dx, c.over_dy, c.over_dz
                    )
                )
                , static_cast<hpx::future<void> >(hpx::when_all(compute_res_futures))
            );

        beginFluid = safe_advance(beginFluid, fluid_cells_.end(), fluid_stride);
        endFluid = safe_advance(endFluid, fluid_cells_.end(), fluid_stride);
    }

    hpx::future<triple<double> > local_max_uv =
        hpx::dataflow(
            hpx::util::unwrapped(
                [](std::vector<triple<double> > max_uvs)
                -> triple<double>
                {
                    triple<double> max_uv(0);

                    for (std::size_t i = 0; i < max_uvs.size(); ++i)
                    {
                        max_uv.x = max_uvs[i].x > max_uv.x ? max_uvs[i].x : max_uv.x;
                        max_uv.y = max_uvs[i].y > max_uv.y ? max_uvs[i].y : max_uv.y;
                        max_uv.z = max_uvs[i].z > max_uv.z ? max_uvs[i].z : max_uv.z;
                    }

                    return max_uv;
                }
            )
            , local_max_uvs
        );

    t_ += dt;
    ++step_;

    return local_max_uv;
}

}
}
}
