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

partition_server::partition_server(io::config const& cfg)
:   c(cfg),
    cells_x_(c.cells_x_per_block * c.num_x_blocks + 2),
    cells_y_(c.cells_y_per_block * c.num_y_blocks + 2),
    cells_z_(c.cells_z_per_block * c.num_z_blocks + 2),
    is_left_(c.idx == 0),
    is_right_(c.idx == c.num_partitions_x - 1),
    is_bottom_(c.idy == 0),
    is_top_(c.idy == c.num_partitions_y - 1),
    is_front_(c.idz == 0),
    is_back_(c.idy == c.num_partitions_z - 1)
{

#ifdef WITH_SOR
    if (c.verbose)
        std::cout << "Solver: SOR" << std::endl;
#else
    if (c.verbose)
        std::cout << "Solver: blockwise Jacobi" << std::endl;
#endif

#ifdef WITH_FOR_EACH
    if (c.verbose)
        std::cout << "Parallelization: hpx::parallel::for_each" << std::endl;
    c.cells_x_per_block = cells_x_ - 2;
    c.cells_y_per_block = cells_y_ - 2;
    c.num_x_blocks = 1;
    c.num_y_blocks = 1;
#else
    if (c.verbose)
        std::cout << "Parellelization: custom grain size" << std::endl;
#endif

    if (c.verbose)
        std::cout << c << std::endl;

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
    fluid_cells_.resize(c.num_x_blocks, c.num_y_blocks, c.num_z_blocks);
    boundary_cells_.resize(c.num_x_blocks, c.num_y_blocks, c.num_z_blocks);
    obstacle_cells_.resize(c.num_x_blocks, c.num_y_blocks, c.num_z_blocks);

    for (std::size_t z = 1, nz_block = 0; z < cells_z_ - 1; z += c.cells_z_per_block, ++nz_block)
    {
        range_type z_range(z, std::min(z + c.cells_z_per_block, cells_z_ - 1));

        for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
        {
            range_type y_range(y, std::min(y + c.cells_y_per_block, cells_y_ - 1));

            for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
            {
                range_type x_range(x, std::min(x + c.cells_x_per_block, cells_x_ - 1));

                for (std::size_t k = z_range.first; k < z_range.second; ++k)
                    for (std::size_t j = y_range.first; j < y_range.second; ++j)
                        for (std::size_t i = x_range.first; i < x_range.second; ++i)
                        {
                            cell_type_data_(i, j, k) = std::move(c.flag_grid[k * cells_x_ * cells_y_ + j * cells_x_ + i]);

                          //  data_[P](i, j, k) = cell_type_data_(i, j, k).to_ulong();
                         //   data_[U](i, j, k) = i - 1;
                          //  data_[V](i, j, k) = j - 1;
                          //  data_[W](i, j, k) = k - 1;

                            if (cell_type_data_(i, j, k).test(is_fluid))
                            {
                                fluid_cells_(nx_block, ny_block, nz_block).emplace_back(i, j, k);

                            }
                            else if (cell_type_data_(i, j, k).test(is_boundary))
                                boundary_cells_(nx_block, ny_block, nz_block).emplace_back(i, j, k);
                            else if (!cell_type_data_(i, j, k).none())
                                obstacle_cells_(nx_block, ny_block, nz_block).emplace_back(i, j, k);
                        }

            }
        }
    }

}

void partition_server::init()
{
   // for (std::size_t var = 0; var < NUM_VARIABLES; ++var)
       // data_[var].clear(0);

    rhs_data_.clear(0);

   // step_ = 0;
    t_ = 0;
    next_out_ = -1e-12;
    outcount_ = 0;
    current = 1;
    last = 0;

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

    send_futures_U.resize(NUM_DIRECTIONS);
    send_futures_V.resize(NUM_DIRECTIONS);
    recv_futures_U.resize(NUM_DIRECTIONS);
    recv_futures_V.resize(NUM_DIRECTIONS);
    send_futures_F.resize(NUM_DIRECTIONS);
    send_futures_G.resize(NUM_DIRECTIONS);
    recv_futures_F.resize(NUM_DIRECTIONS);
    recv_futures_G.resize(NUM_DIRECTIONS);
    send_futures_P.resize(NUM_DIRECTIONS);

    set_velocity_futures.resize(c.num_x_blocks, c.num_y_blocks, c.num_z_blocks);
    init_future_grid(recv_futures_U);
    init_future_grid(recv_futures_V);

    compute_fg_futures.resize(c.num_x_blocks, c.num_y_blocks, c.num_z_blocks);
    init_future_grid(recv_futures_F);
    init_future_grid(recv_futures_G);

    compute_rhs_futures.resize(c.num_x_blocks, c.num_y_blocks, c.num_z_blocks);
    compute_res_futures.resize(c.num_x_blocks, c.num_y_blocks, c.num_z_blocks);

    for (auto& a : compute_res_futures)
        a = hpx::make_ready_future(0.);

    set_p_futures.resize(c.num_x_blocks, c.num_y_blocks, c.num_z_blocks);
    recv_futures_P[current].resize(NUM_DIRECTIONS);
    recv_futures_P[last].resize(NUM_DIRECTIONS);
    init_future_grid(recv_futures_P[current]);
    init_future_grid(recv_futures_P[last]);

    for (auto& a : recv_futures_P[last])
        for (auto& b : a)
            b = hpx::make_ready_future();

    sor_cycle_futures[current].resize(c.num_x_blocks, c.num_y_blocks, c.num_z_blocks);
    sor_cycle_futures[last].resize(c.num_x_blocks, c.num_y_blocks, c.num_z_blocks);

    for (auto& a : sor_cycle_futures[last])
        a = hpx::make_ready_future();

    token.reset();
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


triple<Real> partition_server::do_timestep(Real dt)
{
    for (std::size_t nz_block = 0; nz_block < c.num_z_blocks; ++nz_block)
        for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
            for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
            {
                hpx::shared_future<void> calc_future =
                    hpx::async(
                        hpx::util::bind(
                            &stencils<STENCIL_SET_VELOCITY_OBSTACLE>::call,
                            boost::ref(data_[U]), boost::ref(data_[V]), boost::ref(data_[W]),
                            boost::ref(cell_type_data_),
                            boost::ref(boundary_cells_(nx_block, ny_block, nz_block)),
                            boost::ref(obstacle_cells_(nx_block, ny_block, nz_block)),
                            boost::ref(c.bnd_condition)
                        )
                    );

                 set_velocity_futures(nx_block, ny_block, nz_block) = calc_future;
            }

    hpx::wait_all(set_velocity_futures.data_);

    for (std::size_t nz_block = 0; nz_block < c.num_z_blocks; ++nz_block)
        for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
            for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
            {
                hpx::shared_future<void> calc_future =
                    hpx::async(
                        hpx::util::bind(
                            &stencils<STENCIL_COMPUTE_FG>::call,
                            boost::ref(data_[F]), boost::ref(data_[G]), boost::ref(data_[H]),
                            boost::ref(data_[U]), boost::ref(data_[V]), boost::ref(data_[W]),
                            boost::ref(cell_type_data_),
                            boost::ref(boundary_cells_(nx_block, ny_block, nz_block)),
                            boost::ref(obstacle_cells_(nx_block, ny_block, nz_block)),
                            boost::ref(fluid_cells_(nx_block, ny_block, nz_block)),
                            c.re, c.gx, c.gy, c.gz, c.beta, c.dx, c.dy, c.dz,
                            c.dx_sq, c.dy_sq, c.dz_sq, dt, c.alpha
                        )
                    );

                 compute_fg_futures(nx_block, ny_block, nz_block) = calc_future;
            }

    hpx::wait_all(compute_fg_futures.data_);

    for (std::size_t nz_block = 0; nz_block < c.num_z_blocks; ++nz_block)
        for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
            for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
            {
                hpx::shared_future<void> calc_future =
                    hpx::async(
                        hpx::util::bind(
                            &stencils<STENCIL_COMPUTE_RHS>::call,
                            boost::ref(rhs_data_),
                            boost::ref(data_[F]), boost::ref(data_[G]), boost::ref(data_[H]),
                            boost::ref(cell_type_data_),
                            boost::ref(fluid_cells_(nx_block, ny_block, nz_block)),
                            c.dx, c.dy, c.dz, dt
                        )
                    );

                 compute_rhs_futures(nx_block, ny_block, nz_block) = calc_future;
            }

    hpx::wait_all(compute_rhs_futures.data_);

    Real res = c.eps_sq + 1;

    for (std::size_t iter; iter < c.iter_max && res > c.eps_sq; ++iter)
    {
         for (std::size_t nz_block = 0; nz_block < c.num_z_blocks; ++nz_block)
                for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
                    for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
                    {
                        hpx::shared_future<void> calc_future =
                            hpx::async(
                                hpx::util::bind(
                                    &stencils<STENCIL_SET_P_OBSTACLE>::call,
                                    boost::ref(data_[P]),
                                    boost::ref(cell_type_data_),
                                    boost::ref(boundary_cells_(nx_block, ny_block, nz_block)),
                                    boost::ref(obstacle_cells_(nx_block, ny_block, nz_block)),
                                    token
                                )
                            );

                         set_p_futures(nx_block, ny_block, nz_block) = calc_future;
                    }

            hpx::wait_all(set_p_futures.data_);

            for (std::size_t nz_block = 0; nz_block < c.num_z_blocks; ++nz_block)
                for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
                    for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
                    {
                        hpx::shared_future<void> calc_future =
                            hpx::async(
                                hpx::util::bind(
                                    &stencils<STENCIL_JACOBI>::call,
                                    boost::ref(data_[P]),
                                    boost::ref(rhs_data_),
                                    boost::ref(fluid_cells_(nx_block, ny_block, nz_block)),
                                    c.dx_sq, c.dy_sq, c.dz_sq, token
                                )
                            );

                     /*       hpx::shared_future<void> calc_future =
                            hpx::async(
                                hpx::util::bind(
                                    &stencils<STENCIL_SOR>::call,
                                    boost::ref(data_[P]),
                                    boost::ref(rhs_data_),
                                    boost::ref(fluid_cells_(nx_block, ny_block, nz_block)),
                                    c.part1, c.part2, c.dx_sq, c.dy_sq, c.dz_sq, token
                                )
                            );*/

                         sor_cycle_futures[current](nx_block, ny_block, nz_block) = calc_future;
                    }

            hpx::wait_all(sor_cycle_futures[current].data_);

            for (std::size_t nz_block = 0; nz_block < c.num_z_blocks; ++nz_block)
                for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
                    for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
                    {
                        hpx::shared_future<Real> calc_future =
                            hpx::async(
                                hpx::util::bind(
                                    &stencils<STENCIL_COMPUTE_RESIDUAL>::call,
                                    boost::ref(data_[P]),
                                    boost::ref(rhs_data_),
                                    boost::ref(fluid_cells_(nx_block, ny_block, nz_block)),
                                    c.dx_sq, c.dy_sq, c.dz_sq, token
                                )
                            );

                         compute_res_futures(nx_block, ny_block, nz_block) = calc_future;
                    }

            hpx::wait_all(compute_res_futures.data_);

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

            local_residual.wait();
            res = local_residual.get();

            if (res < c.eps_sq || iter == c.iter_max - 1)
                std::cout << "Iter: " << iter + 1 << "\tResidual: " << res << "\tdt: " << dt << std::endl;
    }

    local_max_uvs.clear();
    local_max_uvs.reserve(c.num_x_blocks * c.num_y_blocks * c.num_z_blocks);

    for (std::size_t nz_block = 0; nz_block < c.num_z_blocks; ++nz_block)
        for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
            for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
            {
                    local_max_uvs.push_back(hpx::async(
                        hpx::util::bind(
                            &stencils<STENCIL_UPDATE_VELOCITY>::call,
                            boost::ref(data_[U]), boost::ref(data_[V]), boost::ref(data_[W]),
                            boost::ref(data_[F]), boost::ref(data_[G]), boost::ref(data_[H]),
                            boost::ref(data_[P]),
                            boost::ref(cell_type_data_),
                            boost::ref(fluid_cells_(nx_block, ny_block, nz_block)),
                            dt, c.over_dx, c.over_dy, c.over_dz
                        )
                    ));

            }

    hpx::future<triple<Real> > local_max_uv =
                hpx::dataflow(
                    hpx::util::unwrapped(
                        [](std::vector<triple<Real> > max_uvs)
                        -> triple<Real>
                        {
                            triple<Real> max_uv(0);

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

    local_max_uv.wait();

    io::writer::write_vtk(data_[P], data_[U], data_[V], data_[W], data_[F], data_[G], data_[H], cell_type_data_, c.num_partitions_x, c.num_partitions_y, c.num_partitions_z,
                          c.i_max, c.j_max, c.k_max, c.dx, c.dy, c.dz, outcount_++, c.rank, c.idx, c.idy, c.idz);

   // hpx::wait_all(compute_res_futures.data_);

    return local_max_uv.get();
}

}
}
}
