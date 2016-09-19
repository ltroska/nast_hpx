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
    is_bottom_(c.idz == 0),
    is_top_(c.idz == c.num_partitions_z - 1),
    is_front_(c.idy == 0),
    is_back_(c.idy == c.num_partitions_y - 1)
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

                            if (cell_type_data_(i, j, k).test(is_fluid))
                                fluid_cells_(nx_block, ny_block, nz_block).emplace_back(i, j, k);
                            else if (cell_type_data_(i, j, k).count() > 2)
                                obstacle_cells_(nx_block, ny_block, nz_block).emplace_back(i, j, k);
                        }

            }
        }
    }

    std::cout << boundary_cells_.data_.size() << std::endl;
}

void partition_server::init()
{
    for (std::size_t var = 0; var < NUM_VARIABLES; ++var)
        data_[var].clear(0);

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
        send_buffer_left_.dest_ = ids_[c.idz * c.num_partitions_x * c.num_partitions_y + c.idy * c.num_partitions_x + c.idx - 1];
        recv_buffer_left_[U].valid_ = true;
        recv_buffer_left_[V].valid_ = true;
        recv_buffer_left_[W].valid_ = true;
        recv_buffer_left_[F].valid_ = true;
        recv_buffer_left_[P].valid_ = true;
    }

    if (!is_right_)
    {
        send_buffer_right_.dest_ = ids_[c.idz * c.num_partitions_x * c.num_partitions_y + c.idy * c.num_partitions_x + c.idx + 1];

        recv_buffer_right_[U].valid_ = true;
        recv_buffer_right_[V].valid_ = true;
        recv_buffer_right_[W].valid_ = true;
        recv_buffer_right_[P].valid_ = true;
    }

    if (!is_bottom_)
    {
        send_buffer_bottom_.dest_ = ids_[(c.idz - 1) * c.num_partitions_x * c.num_partitions_y + c.idy * c.num_partitions_x + c.idx];

        recv_buffer_right_[U].valid_ = true;
        recv_buffer_right_[V].valid_ = true;
        recv_buffer_right_[W].valid_ = true;
        recv_buffer_right_[H].valid_ = true;
        recv_buffer_right_[P].valid_ = true;    }

    if (!is_top_)
    {
        send_buffer_top_.dest_ = ids_[(c.idz + 1) * c.num_partitions_x * c.num_partitions_y + c.idy * c.num_partitions_x + c.idx];

        recv_buffer_right_[U].valid_ = true;
        recv_buffer_right_[V].valid_ = true;
        recv_buffer_right_[W].valid_ = true;
        recv_buffer_right_[P].valid_ = true;    }

    if (!is_front_)
    {
        send_buffer_front_.dest_ = ids_[c.idz * c.num_partitions_x * c.num_partitions_y + (c.idy - 1) * c.num_partitions_x + c.idx];

        recv_buffer_right_[U].valid_ = true;
        recv_buffer_right_[V].valid_ = true;
        recv_buffer_right_[W].valid_ = true;
        recv_buffer_right_[G].valid_ = true;
        recv_buffer_right_[P].valid_ = true;    }

    if (!is_back_)
    {
        send_buffer_back_.dest_ = ids_[c.idz * c.num_partitions_x * c.num_partitions_y + (c.idy + 1) * c.num_partitions_x + c.idx];

        recv_buffer_right_[U].valid_ = true;
        recv_buffer_right_[V].valid_ = true;
        recv_buffer_right_[W].valid_ = true;
        recv_buffer_right_[P].valid_ = true;    }

    if (!is_back_ && !is_left_)
    {
        send_buffer_back_left_.dest_ = ids_[c.idz * c.num_partitions_x * c.num_partitions_y + (c.idy + 1) * c.num_partitions_x + c.idx - 1];

        recv_buffer_back_left_[U].valid_ = true;
    }

    if (!is_front_ && !is_right_)
    {
        send_buffer_front_right_.dest_ = ids_[c.idz * c.num_partitions_x * c.num_partitions_y + (c.idy - 1) * c.num_partitions_x + c.idx + 1];

        recv_buffer_front_right_[V].valid_ = true;
    }

    if (!is_bottom_ && !is_right_)
    {
        send_buffer_bottom_right_.dest_ = ids_[(c.idz - 1) * c.num_partitions_x * c.num_partitions_y + c.idy * c.num_partitions_x + c.idx + 1];

        recv_buffer_bottom_right_[W].valid_ = true;
    }

    if (!is_top_ && !is_left_)
    {
        send_buffer_top_left_.dest_ = ids_[(c.idz + 1) * c.num_partitions_x * c.num_partitions_y + c.idy * c.num_partitions_x + c.idx - 1];

        recv_buffer_top_left_[U].valid_ = true;
    }

    if (!is_back_ && !is_bottom_)
    {
        send_buffer_back_bottom_.dest_ = ids_[(c.idz - 1) * c.num_partitions_x * c.num_partitions_y + (c.idy + 1) * c.num_partitions_x + c.idx];

        recv_buffer_back_bottom_[W].valid_ = true;
    }

    if (!is_front_ && !is_top_)
    {
        send_buffer_front_top_.dest_ = ids_[(c.idz + 1) * c.num_partitions_x * c.num_partitions_y + (c.idy - 1) * c.num_partitions_x + c.idx];

        recv_buffer_front_top_[V].valid_ = true;
    }


    recv_futures.resize(NUM_VARIABLES);

    for (std::size_t var = 0; var < NUM_VARIABLES; ++var)
    {
        recv_futures[var].resize(NUM_DIRECTIONS);

        recv_futures[var][RIGHT].resize(c.num_y_blocks, c.num_z_blocks);
        recv_futures[var][LEFT].resize(c.num_y_blocks, c.num_z_blocks);
        recv_futures[var][BOTTOM].resize(c.num_x_blocks, c.num_y_blocks);
        recv_futures[var][TOP].resize(c.num_x_blocks, c.num_y_blocks);
        recv_futures[var][FRONT].resize(c.num_x_blocks, c.num_z_blocks);
        recv_futures[var][BACK].resize(c.num_x_blocks, c.num_z_blocks);
        recv_futures[var][BACK_LEFT].resize(c.num_z_blocks, 1);
        recv_futures[var][FRONT_RIGHT].resize(c.num_z_blocks, 1);
        recv_futures[var][BOTTOM_RIGHT].resize(c.num_y_blocks, 1);
        recv_futures[var][TOP_LEFT].resize(c.num_y_blocks, 1);
        recv_futures[var][BACK_BOTTOM].resize(c.num_x_blocks, 1);
        recv_futures[var][FRONT_TOP].resize(c.num_x_blocks, 1);
    }

    set_velocity_futures.resize(c.num_x_blocks, c.num_y_blocks, c.num_z_blocks);

    compute_fg_futures.resize(c.num_x_blocks, c.num_y_blocks, c.num_z_blocks);


    compute_rhs_futures.resize(c.num_x_blocks, c.num_y_blocks, c.num_z_blocks);
    compute_res_futures.resize(c.num_x_blocks, c.num_y_blocks, c.num_z_blocks);

    for (auto& a : compute_res_futures)
        a = hpx::make_ready_future(0.);

    set_p_futures.resize(c.num_x_blocks, c.num_y_blocks, c.num_z_blocks);

    solver_cycle_futures.resize(c.num_x_blocks, c.num_y_blocks, c.num_z_blocks);

    for (auto& a : solver_cycle_futures)
        a = hpx::make_ready_future();

    token.reset();
}

template<>
void partition_server::send_boundary<LEFT>(std::size_t step, std::size_t var, partition_data<hpx::shared_future<void> >& send_futures)
{
    for (std::size_t y = 0; y < send_futures.size_y_; ++y)
        for (std::size_t z = 0; z < send_futures.size_z_; ++z)
            send_futures(0, y, z).then(
                hpx::launch::async,
                hpx::util::bind(
                    boost::ref(send_buffer_left_),
                    boost::ref(data_[var]),
                    step * c.num_y_blocks * c.num_z_blocks + (z * c.num_y_blocks + y),
                    var,
                    y,
                    z,
                    c.cells_y_per_block,
                    c.cells_z_per_block
                )
            );
}

template<>
void partition_server::send_boundary<RIGHT>(std::size_t step, std::size_t var, partition_data<hpx::shared_future<void> >& send_futures)
{
    for (std::size_t y = 0; y < send_futures.size_y_; ++y)
        for (std::size_t z = 0; z < send_futures.size_z_; ++z)
            send_futures(send_futures.size_x_ - 1, y, z).then(
                hpx::launch::async,
                hpx::util::bind(
                    boost::ref(send_buffer_right_),
                    boost::ref(data_[var]),
                    step * c.num_y_blocks * c.num_z_blocks + (z * c.num_y_blocks + y),
                    var,
                    y,
                    z,
                    c.cells_y_per_block,
                    c.cells_z_per_block
                )
            );
}

template<>
void partition_server::send_boundary<BOTTOM>(std::size_t step, std::size_t var, partition_data<hpx::shared_future<void> >& send_futures)
{
    for (std::size_t x = 0; x < send_futures.size_x_; ++x)
        for (std::size_t y = 0; y < send_futures.size_y_; ++y)
            send_futures(x, y, 0).then(
                hpx::launch::async,
                hpx::util::bind(
                    boost::ref(send_buffer_bottom_),
                    boost::ref(data_[var]),
                    step * c.num_x_blocks * c.num_y_blocks + (y * c.num_x_blocks + x),
                    var,
                    x,
                    y,
                    c.cells_x_per_block,
                    c.cells_y_per_block
                )
            );
}

template<>
void partition_server::send_boundary<TOP>(std::size_t step, std::size_t var, partition_data<hpx::shared_future<void> >& send_futures)
{
    for (std::size_t x = 0; x < send_futures.size_x_; ++x)
        for (std::size_t y = 0; y < send_futures.size_y_; ++y)
            send_futures(x, y, send_futures.size_z_ - 1).then(
                hpx::launch::async,
                hpx::util::bind(
                    boost::ref(send_buffer_top_),
                    boost::ref(data_[var]),
                    step * c.num_x_blocks * c.num_y_blocks + (y * c.num_x_blocks + x),
                    var,
                    x,
                    y,
                    c.cells_x_per_block,
                    c.cells_y_per_block
                )
            );
}

template<>
void partition_server::send_boundary<FRONT>(std::size_t step, std::size_t var, partition_data<hpx::shared_future<void> >& send_futures)
{
    for (std::size_t x = 0; x < send_futures.size_x_; ++x)
        for (std::size_t z = 0; z < send_futures.size_z_; ++z)
            send_futures(x, 0, z).then(
                hpx::launch::async,
                hpx::util::bind(
                    boost::ref(send_buffer_front_),
                    boost::ref(data_[var]),
                    step * c.num_x_blocks * c.num_z_blocks + (z * c.num_x_blocks + x),
                    var,
                    x,
                    z,
                    c.cells_x_per_block,
                    c.cells_z_per_block
                )
            );
}

template<>
void partition_server::send_boundary<BACK>(std::size_t step, std::size_t var, partition_data<hpx::shared_future<void> >& send_futures)
{
    for (std::size_t x = 0; x < send_futures.size_x_; ++x)
        for (std::size_t z = 0; z < send_futures.size_z_; ++z)
            send_futures(x, send_futures.size_y_ - 1, z).then(
                hpx::launch::async,
                hpx::util::bind(
                    boost::ref(send_buffer_back_),
                    boost::ref(data_[var]),
                    step * c.num_x_blocks * c.num_z_blocks + (z * c.num_x_blocks + x),
                    var,
                    x,
                    z,
                    c.cells_x_per_block,
                    c.cells_z_per_block
                )
            );
}

template<>
void partition_server::send_boundary<BACK_LEFT>(std::size_t step, std::size_t var, partition_data<hpx::shared_future<void> >& send_futures)
{
    for (std::size_t z = 0; z < send_futures.size_z_; ++z)
        send_futures(0, send_futures.size_y_ - 1, z).then(
            hpx::launch::async,
            hpx::util::bind(
                boost::ref(send_buffer_back_left_),
                boost::ref(data_[var]),
                step * c.num_z_blocks + z,
                var,
                z,
                0,
                c.cells_z_per_block,
                0
            )
        );
}

template<>
void partition_server::send_boundary<FRONT_RIGHT>(std::size_t step, std::size_t var, partition_data<hpx::shared_future<void> >& send_futures)
{
    for (std::size_t z = 0; z < send_futures.size_z_; ++z)
        send_futures(send_futures.size_x_ - 1, 0, z).then(
            hpx::launch::async,
            hpx::util::bind(
                boost::ref(send_buffer_front_right_),
                boost::ref(data_[var]),
                step * c.num_z_blocks + z,
                var,
                z,
                0,
                c.cells_z_per_block,
                0
            )
        );
}

template<>
void partition_server::send_boundary<BOTTOM_RIGHT>(std::size_t step, std::size_t var, partition_data<hpx::shared_future<void> >& send_futures)
{
    for (std::size_t y = 0; y < send_futures.size_y_; ++y)
        send_futures(send_futures.size_x_ - 1, y, 0).then(
            hpx::launch::async,
            hpx::util::bind(
                boost::ref(send_buffer_bottom_right_),
                boost::ref(data_[var]),
                step * c.num_y_blocks + y,
                var,
                y,
                0,
                c.cells_y_per_block,
                0
            )
        );
}

template<>
void partition_server::send_boundary<TOP_LEFT>(std::size_t step, std::size_t var, partition_data<hpx::shared_future<void> >& send_futures)
{
    for (std::size_t y = 0; y < send_futures.size_y_; ++y)
        send_futures(0, y, send_futures.size_z_ - 1).then(
            hpx::launch::async,
            hpx::util::bind(
                boost::ref(send_buffer_top_left_),
                boost::ref(data_[var]),
                step * c.num_y_blocks + y,
                var,
                y,
                0,
                c.cells_y_per_block,
                0
            )
        );
}

template<>
void partition_server::send_boundary<BACK_BOTTOM>(std::size_t step, std::size_t var, partition_data<hpx::shared_future<void> >& send_futures)
{
    for (std::size_t x = 0; x < send_futures.size_x_; ++x)
        send_futures(x, send_futures.size_y_ - 1, 0).then(
            hpx::launch::async,
            hpx::util::bind(
                boost::ref(send_buffer_back_bottom_),
                boost::ref(data_[var]),
                step * c.num_x_blocks + x,
                var,
                x,
                0,
                c.cells_x_per_block,
                0
            )
        );
}

template<>
void partition_server::send_boundary<FRONT_TOP>(std::size_t step, std::size_t var, partition_data<hpx::shared_future<void> >& send_futures)
{
    for (std::size_t x = 0; x < send_futures.size_x_; ++x)
        send_futures(x, 0, send_futures.size_z_ - 1).then(
            hpx::launch::async,
            hpx::util::bind(
                boost::ref(send_buffer_front_top_),
                boost::ref(data_[var]),
                step * c.num_x_blocks + x,
                var,
                x,
                0,
                c.cells_x_per_block,
                0
            )
        );
}

template<>
void partition_server::receive_boundary<LEFT>(std::size_t step, std::size_t var, future_grid_4d& recv_futures)
{
     for (std::size_t y = 0; y < c.num_y_blocks; ++y)
        for (std::size_t z = 0; z < c.num_z_blocks; ++z)
            recv_futures[var][LEFT](y, z) =
                hpx::async(
                    hpx::util::bind(
                        boost::ref(recv_buffer_left_[var]),
                        boost::ref(data_[var]),
                        step * c.num_y_blocks * c.num_z_blocks + (z * c.num_y_blocks + y),
                        var,
                        y,
                        z,
                        c.cells_y_per_block,
                        c.cells_z_per_block
                    )
                );
}

template<>
void partition_server::receive_boundary<RIGHT>(std::size_t step, std::size_t var, future_grid_4d& recv_futures)
{
     for (std::size_t y = 0; y < c.num_y_blocks; ++y)
        for (std::size_t z = 0; z < c.num_z_blocks; ++z)
            recv_futures[var][RIGHT](y, z) =
                hpx::async(
                    hpx::util::bind(
                        boost::ref(recv_buffer_right_[var]),
                        boost::ref(data_[var]),
                        step * c.num_y_blocks * c.num_z_blocks + (z * c.num_y_blocks + y),
                        var,
                        y,
                        z,
                        c.cells_y_per_block,
                        c.cells_z_per_block
                    )
                );
}

template<>
void partition_server::receive_boundary<BOTTOM>(std::size_t step, std::size_t var, future_grid_4d& recv_futures)
{
     for (std::size_t x = 0; x < c.num_x_blocks; ++x)
        for (std::size_t y = 0; y < c.num_y_blocks; ++y)
            recv_futures[var][BOTTOM](x, y) =
                hpx::async(
                    hpx::util::bind(
                        boost::ref(recv_buffer_bottom_[var]),
                        boost::ref(data_[var]),
                        step * c.num_x_blocks * c.num_y_blocks + (y * c.num_x_blocks + x),
                        var,
                        x,
                        y,
                        c.cells_x_per_block,
                        c.cells_y_per_block
                    )
                );
}

template<>
void partition_server::receive_boundary<TOP>(std::size_t step, std::size_t var, future_grid_4d& recv_futures)
{
     for (std::size_t x = 0; x < c.num_x_blocks; ++x)
        for (std::size_t y = 0; y < c.num_y_blocks; ++y)
            recv_futures[var][TOP](x, y) =
                hpx::async(
                    hpx::util::bind(
                        boost::ref(recv_buffer_top_[var]),
                        boost::ref(data_[var]),
                        step * c.num_x_blocks * c.num_y_blocks + (y * c.num_x_blocks + x),
                        var,
                        x,
                        y,
                        c.cells_x_per_block,
                        c.cells_y_per_block
                    )
                );
}

template<>
void partition_server::receive_boundary<FRONT>(std::size_t step, std::size_t var, future_grid_4d& recv_futures)
{
     for (std::size_t x = 0; x < c.num_x_blocks; ++x)
        for (std::size_t z = 0; z < c.num_z_blocks; ++z)
            recv_futures[var][FRONT](x, z) =
                hpx::async(
                    hpx::util::bind(
                        boost::ref(recv_buffer_front_[var]),
                        boost::ref(data_[var]),
                        step * c.num_x_blocks * c.num_z_blocks + (z * c.num_x_blocks + x),
                        var,
                        x,
                        z,
                        c.cells_x_per_block,
                        c.cells_z_per_block
                    )
                );
}

template<>
void partition_server::receive_boundary<BACK>(std::size_t step, std::size_t var, future_grid_4d& recv_futures)
{
     for (std::size_t x = 0; x < c.num_x_blocks; ++x)
        for (std::size_t z = 0; z < c.num_z_blocks; ++z)
            recv_futures[var][BACK](x, z) =
                hpx::async(
                    hpx::util::bind(
                        boost::ref(recv_buffer_back_[var]),
                        boost::ref(data_[var]),
                        step * c.num_x_blocks * c.num_z_blocks + (z * c.num_x_blocks + x),
                        var,
                        x,
                        z,
                        c.cells_x_per_block,
                        c.cells_z_per_block
                    )
                );
}

template<>
void partition_server::receive_boundary<BACK_LEFT>(std::size_t step, std::size_t var, future_grid_4d& recv_futures)
{
    for (std::size_t z = 0; z < c.num_z_blocks; ++z)
        recv_futures[var][BACK_LEFT](z, 0) =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_back_left_[var]),
                    boost::ref(data_[var]),
                    step * c.num_z_blocks + z,
                    var,
                    z,
                    0,
                    c.cells_z_per_block,
                    0
                )
            );
}

template<>
void partition_server::receive_boundary<FRONT_RIGHT>(std::size_t step, std::size_t var, future_grid_4d& recv_futures)
{
    for (std::size_t z = 0; z < c.num_z_blocks; ++z)
        recv_futures[var][FRONT_RIGHT](z, 0) =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_front_right_[var]),
                    boost::ref(data_[var]),
                    step * c.num_z_blocks + z,
                    var,
                    z,
                    0,
                    c.cells_z_per_block,
                    0
                )
            );
}

template<>
void partition_server::receive_boundary<BOTTOM_RIGHT>(std::size_t step, std::size_t var, future_grid_4d& recv_futures)
{
    for (std::size_t y = 0; y < c.num_y_blocks; ++y)
        recv_futures[var][BOTTOM_RIGHT](y, 0) =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_bottom_right_[var]),
                    boost::ref(data_[var]),
                    step * c.num_y_blocks + y,
                    var,
                    y,
                    0,
                    c.cells_y_per_block,
                    0
                )
            );
}

template<>
void partition_server::receive_boundary<TOP_LEFT>(std::size_t step, std::size_t var, future_grid_4d& recv_futures)
{
    for (std::size_t y = 0; y < c.num_y_blocks; ++y)
        recv_futures[var][TOP_LEFT](y, 0) =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_top_left_[var]),
                    boost::ref(data_[var]),
                    step * c.num_y_blocks + y,
                    var,
                    y,
                    0,
                    c.cells_y_per_block,
                    0
                )
            );
}

template<>
void partition_server::receive_boundary<BACK_BOTTOM>(std::size_t step, std::size_t var, future_grid_4d& recv_futures)
{
    for (std::size_t x = 0; x < c.num_x_blocks; ++x)
        recv_futures[var][BACK_BOTTOM](x, 0) =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_back_bottom_[var]),
                    boost::ref(data_[var]),
                    step * c.num_x_blocks + x,
                    var,
                    x,
                    0,
                    c.cells_x_per_block,
                    0
                )
            );
}

template<>
void partition_server::receive_boundary<FRONT_TOP>(std::size_t step, std::size_t var, future_grid_4d& recv_futures)
{
    for (std::size_t x = 0; x < c.num_x_blocks; ++x)
        recv_futures[var][FRONT_TOP](x, 0) =
            hpx::async(
                hpx::util::bind(
                    boost::ref(recv_buffer_front_top_[var]),
                    boost::ref(data_[var]),
                    step * c.num_x_blocks + x,
                    var,
                    x,
                    0,
                    c.cells_x_per_block,
                    0
                )
            );
}

template<>
hpx::shared_future<void> partition_server::get_dependency<LEFT>(std::size_t idx_block, std::size_t idy_block, std::size_t idz_block,
    future_grid_3d const& recv_futures, partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if (is_left_ && idx_block == 0)
        return hpx::make_ready_future();
    else if (idx_block == 0)
        return recv_futures[LEFT](idy_block, idz_block);

    return calc_futures(idx_block - 1, idy_block, idz_block);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<RIGHT>(std::size_t idx_block, std::size_t idy_block, std::size_t idz_block,
    future_grid_3d const& recv_futures, partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if (is_right_ && idx_block == calc_futures.size_x_ - 1)
        return hpx::make_ready_future();
    else if (idx_block == calc_futures.size_x_ - 1)
        return recv_futures[RIGHT](idy_block, idz_block);

    return calc_futures(idx_block + 1, idy_block, idz_block);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<BOTTOM>(std::size_t idx_block, std::size_t idy_block, std::size_t idz_block,
    future_grid_3d const& recv_futures, partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if (is_bottom_ && idz_block == 0)
        return hpx::make_ready_future();
    else if (idz_block == 0)
        return recv_futures[BOTTOM](idx_block, idy_block);

    return calc_futures(idx_block, idy_block, idz_block - 1);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<TOP>(std::size_t idx_block, std::size_t idy_block, std::size_t idz_block,
    future_grid_3d const& recv_futures, partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if (is_top_ && idz_block == calc_futures.size_z_ - 1)
        return hpx::make_ready_future();
    else if (idz_block == calc_futures.size_z_ - 1)
        return recv_futures[TOP](idx_block, idy_block);

    return calc_futures(idx_block, idy_block, idz_block + 1);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<FRONT>(std::size_t idx_block, std::size_t idy_block, std::size_t idz_block,
    future_grid_3d const& recv_futures, partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if (is_front_ && idy_block == 0)
        return hpx::make_ready_future();
    else if (idy_block == 0)
        return recv_futures[FRONT](idx_block, idz_block);

    return calc_futures(idx_block, idy_block - 1, idz_block);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<BACK>(std::size_t idx_block, std::size_t idy_block, std::size_t idz_block,
    future_grid_3d const& recv_futures, partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if (is_back_ && idy_block == calc_futures.size_y_ - 1)
        return hpx::make_ready_future();
    else if (idy_block == calc_futures.size_y_ - 1)
        return recv_futures[BACK](idx_block, idz_block);

    return calc_futures(idx_block, idy_block + 1, idz_block);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<BACK_LEFT>(std::size_t idx_block, std::size_t idy_block, std::size_t idz_block,
    future_grid_3d const& recv_futures, partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if ((is_back_ && idy_block == calc_futures.size_y_ - 1) || (is_left_ && idx_block == 0))
        return hpx::make_ready_future();

    if (idx_block == 0 && idy_block == calc_futures.size_y_ - 1)
        return recv_futures[BACK_LEFT](idz_block, 0);
    else if (idx_block == 0)
        return recv_futures[LEFT](idy_block + 1, idz_block);
    else if (idy_block == calc_futures.size_y_ - 1)
        return recv_futures[BACK](idx_block - 1, idz_block);

    return calc_futures(idx_block - 1, idy_block + 1, idz_block);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<FRONT_RIGHT>(std::size_t idx_block, std::size_t idy_block, std::size_t idz_block,
    future_grid_3d const& recv_futures, partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if ((is_front_ && idy_block == 0) || (is_right_ && idx_block == calc_futures.size_x_ - 1))
        return hpx::make_ready_future();

    if (idx_block == calc_futures.size_x_ - 1 && idy_block == 0)
        return recv_futures[FRONT_RIGHT](idz_block, 0);
    else if (idx_block == calc_futures.size_x_ - 1)
        return recv_futures[RIGHT](idy_block - 1, idz_block);
    else if (idy_block == 0)
        return recv_futures[FRONT](idx_block + 1, idz_block);

    return calc_futures(idx_block + 1, idy_block - 1, idz_block);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<BOTTOM_RIGHT>(std::size_t idx_block, std::size_t idy_block, std::size_t idz_block,
    future_grid_3d const& recv_futures, partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if ((is_bottom_ && idz_block == 0) || (is_right_ && idx_block == calc_futures.size_x_ - 1))
        return hpx::make_ready_future();

    if (idz_block == 0 && idx_block == calc_futures.size_x_ - 1)
        return recv_futures[BOTTOM_RIGHT](idy_block, 0);
    else if (idz_block == 0)
        return recv_futures[BOTTOM](idx_block + 1, idy_block);
    else if (idx_block == calc_futures.size_x_ - 1)
        return recv_futures[RIGHT](idy_block, idz_block - 1);

    return calc_futures(idx_block + 1, idy_block, idz_block - 1);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<TOP_LEFT>(std::size_t idx_block, std::size_t idy_block, std::size_t idz_block,
    future_grid_3d const& recv_futures, partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if ((is_top_ && idz_block == calc_futures.size_z_ - 1) || (is_left_ && idx_block == 0))
        return hpx::make_ready_future();

    if (idx_block == 0 && idz_block == calc_futures.size_z_ - 1)
        return recv_futures[TOP_LEFT](idy_block, 0);
    else if (idx_block == 0)
        return recv_futures[LEFT](idy_block, idz_block + 1);
    else if (idz_block == calc_futures.size_z_ - 1)
        return recv_futures[TOP](idx_block - 1, idy_block);

    return calc_futures(idx_block - 1, idy_block, idz_block + 1);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<BACK_BOTTOM>(std::size_t idx_block, std::size_t idy_block, std::size_t idz_block,
    future_grid_3d const& recv_futures, partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if ((is_back_ && idy_block == calc_futures.size_y_ - 1) || (is_bottom_ && idz_block == 0))
        return hpx::make_ready_future();

    if (idy_block == calc_futures.size_y_ - 1 && idz_block == 0)
        return recv_futures[BACK_BOTTOM](idx_block, 0);
    else if (idz_block == 0)
        return recv_futures[BOTTOM](idx_block, idy_block + 1);
    else if (idy_block == calc_futures.size_y_ - 1)
        return recv_futures[BACK](idx_block, idz_block - 1);

    return calc_futures(idx_block, idy_block + 1, idz_block - 1);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<FRONT_TOP>(std::size_t idx_block, std::size_t idy_block, std::size_t idz_block,
    future_grid_3d const& recv_futures, partition_data<hpx::shared_future<void> > const& calc_futures)
{
    if ((is_front_ && idy_block == 0) || (is_top_ && idz_block == calc_futures.size_z_ - 1))
        return hpx::make_ready_future();

    if (idy_block == 0 && idz_block == calc_futures.size_z_ - 1)
        return recv_futures[FRONT_TOP](idx_block, 0);
    else if (idy_block == 0)
        return recv_futures[FRONT](idx_block, idz_block + 1);
    else if (idz_block == calc_futures.size_z_ - 1)
        return recv_futures[TOP](idx_block, idy_block - 1);

    return calc_futures(idx_block, idy_block - 1, idz_block + 1);
}

void partition_server::send_boundaries_U(partition_data<hpx::shared_future<void> >& send_futures, std::size_t step)
{
    if (!is_left_)
        send_boundary<LEFT>(step, U, send_futures);

    if (!is_right_)
        send_boundary<RIGHT>(step, U, send_futures);

    if (!is_bottom_)
        send_boundary<BOTTOM>(step, U, send_futures);

    if (!is_top_)
        send_boundary<TOP>(step, U, send_futures);

    if (!is_front_)
        send_boundary<FRONT>(step, U, send_futures);

    if (!is_back_)
        send_boundary<BACK>(step, U, send_futures);

    if (!is_front_ && !is_right_)
        send_boundary<FRONT_RIGHT>(step, U, send_futures);

    if (!is_bottom_ && !is_right_)
        send_boundary<BOTTOM_RIGHT>(step, U, send_futures);
}

void partition_server::receive_boundaries_U(future_grid_4d& recv_futures, std::size_t step)
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

void partition_server::send_boundaries_V(partition_data<hpx::shared_future<void> >& send_futures, std::size_t step)
{
    if (!is_left_)
        send_boundary<LEFT>(step, V, send_futures);

    if (!is_right_)
        send_boundary<RIGHT>(step, V, send_futures);

    if (!is_bottom_)
        send_boundary<BOTTOM>(step, V, send_futures);

    if (!is_top_)
        send_boundary<TOP>(step, V, send_futures);

    if (!is_front_)
        send_boundary<FRONT>(step, V, send_futures);

    if (!is_back_)
        send_boundary<BACK>(step, V, send_futures);

    if (!is_back_ && !is_left_)
        send_boundary<BACK_LEFT>(step, V, send_futures);

    if (!is_back_ && !is_bottom_)
        send_boundary<BACK_BOTTOM>(step, V, send_futures);
}

void partition_server::receive_boundaries_V(future_grid_4d& recv_futures, std::size_t step)
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

void partition_server::send_boundaries_W(partition_data<hpx::shared_future<void> >& send_futures, std::size_t step)
{
    if (!is_left_)
        send_boundary<LEFT>(step, W, send_futures);

    if (!is_right_)
        send_boundary<RIGHT>(step, W, send_futures);

    if (!is_bottom_)
        send_boundary<BOTTOM>(step, W, send_futures);

    if (!is_top_)
        send_boundary<TOP>(step, W, send_futures);

    if (!is_front_)
        send_boundary<FRONT>(step, W, send_futures);

    if (!is_back_)
        send_boundary<BACK>(step, W, send_futures);

    if (!is_top_ && !is_left_)
        send_boundary<TOP_LEFT>(step, W, send_futures);

    if (!is_top_ && !is_front_)
        send_boundary<FRONT_TOP>(step, W, send_futures);
}

void partition_server::receive_boundaries_W(future_grid_4d& recv_futures, std::size_t step)
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

void partition_server::send_boundaries_F(partition_data<hpx::shared_future<void> >& send_futures, std::size_t step)
{
    if (!is_right_)
        send_boundary<RIGHT>(step, F, send_futures);
}

void partition_server::receive_boundaries_F(future_grid_4d& recv_futures, std::size_t step)
{
    if (!is_left_)
        receive_boundary<LEFT>(step, F, recv_futures);
}

void partition_server::send_boundaries_G(partition_data<hpx::shared_future<void> >& send_futures, std::size_t step)
{
    if (!is_back_)
        send_boundary<BACK>(step, G, send_futures);
}

void partition_server::receive_boundaries_G(future_grid_4d& recv_futures, std::size_t step)
{
    if (!is_front_)
        receive_boundary<FRONT>(step, G, recv_futures);
}

void partition_server::send_boundaries_H(partition_data<hpx::shared_future<void> >& send_futures, std::size_t step)
{
    if (!is_top_)
        send_boundary<TOP>(step, H, send_futures);
}


void partition_server::receive_boundaries_H(future_grid_4d& recv_futures, std::size_t step)
{
    if (!is_bottom_)
        receive_boundary<BOTTOM>(step, H, recv_futures);
}

void partition_server::send_boundaries_P(partition_data<hpx::shared_future<void> >& send_futures, std::size_t step)
{
    if (!is_left_)
        send_boundary<LEFT>(step, P, send_futures);

    if (!is_right_)
        send_boundary<RIGHT>(step, P, send_futures);

    if (!is_bottom_)
        send_boundary<BOTTOM>(step, P, send_futures);

    if (!is_top_)
        send_boundary<TOP>(step, P, send_futures);

    if (!is_front_)
        send_boundary<FRONT>(step, P, send_futures);

    if (!is_back_)
        send_boundary<BACK>(step, P, send_futures);
}

void partition_server::receive_boundaries_P(future_grid_4d& recv_futures, std::size_t step)
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

    send_boundaries_U(set_velocity_futures, step_);
    receive_boundaries_U(recv_futures, step_);

    send_boundaries_V(set_velocity_futures, step_);
    receive_boundaries_V(recv_futures, step_);

    send_boundaries_W(set_velocity_futures, step_);
    receive_boundaries_W(recv_futures, step_);

    if (c.vtk && next_out_ < t_)
    {
        next_out_ += c.delta_vec;

        if (c.verbose)
            std::cout << "Output to .vtk in step " << step_ << std::endl;

        hpx::when_all(set_velocity_futures.data_).then(
            hpx::launch::async,
            hpx::util::bind(
                &io::writer::write_vtk,
                boost::ref(data_[P]), boost::ref(data_[U]), boost::ref(data_[V]), boost::ref(data_[W]), boost::ref(cell_type_data_),
                c.num_partitions_x, c.num_partitions_y, c.num_partitions_z, c.i_max, c.j_max, c.k_max, c.dx, c.dx, c.dz, outcount_++,
                c.rank, c.idx, c.idy, c.idz
            )
        ).wait();
    }

    for (std::size_t nz_block = 0; nz_block < c.num_z_blocks; ++nz_block)
        for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
            for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
            {
                hpx::shared_future<void> calc_future =
                    hpx::dataflow(
                        hpx::util::unwrapped(
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
                        )
                        , set_velocity_futures(nx_block, ny_block, nz_block)
                        , get_dependency<LEFT>(nx_block, ny_block, nz_block, recv_futures[U], set_velocity_futures)
                        , get_dependency<RIGHT>(nx_block, ny_block, nz_block, recv_futures[U], set_velocity_futures)
                        , get_dependency<BOTTOM>(nx_block, ny_block, nz_block, recv_futures[U], set_velocity_futures)
                        , get_dependency<TOP>(nx_block, ny_block, nz_block, recv_futures[U], set_velocity_futures)
                        , get_dependency<FRONT>(nx_block, ny_block, nz_block, recv_futures[U], set_velocity_futures)
                        , get_dependency<BACK>(nx_block, ny_block, nz_block, recv_futures[U], set_velocity_futures)
                        , get_dependency<BACK_LEFT>(nx_block, ny_block, nz_block, recv_futures[U], set_velocity_futures)
                        , get_dependency<TOP_LEFT>(nx_block, ny_block, nz_block, recv_futures[U], set_velocity_futures)
                        , get_dependency<LEFT>(nx_block, ny_block, nz_block, recv_futures[V], set_velocity_futures)
                        , get_dependency<RIGHT>(nx_block, ny_block, nz_block, recv_futures[V], set_velocity_futures)
                        , get_dependency<BOTTOM>(nx_block, ny_block, nz_block, recv_futures[V], set_velocity_futures)
                        , get_dependency<TOP>(nx_block, ny_block, nz_block, recv_futures[V], set_velocity_futures)
                        , get_dependency<FRONT>(nx_block, ny_block, nz_block, recv_futures[V], set_velocity_futures)
                        , get_dependency<BACK>(nx_block, ny_block, nz_block, recv_futures[V], set_velocity_futures)
                        , get_dependency<FRONT_RIGHT>(nx_block, ny_block, nz_block, recv_futures[V], set_velocity_futures)
                        , get_dependency<FRONT_TOP>(nx_block, ny_block, nz_block, recv_futures[V], set_velocity_futures)
                        , get_dependency<LEFT>(nx_block, ny_block, nz_block, recv_futures[W], set_velocity_futures)
                        , get_dependency<RIGHT>(nx_block, ny_block, nz_block, recv_futures[W], set_velocity_futures)
                        , get_dependency<BOTTOM>(nx_block, ny_block, nz_block, recv_futures[W], set_velocity_futures)
                        , get_dependency<TOP>(nx_block, ny_block, nz_block, recv_futures[W], set_velocity_futures)
                        , get_dependency<FRONT>(nx_block, ny_block, nz_block, recv_futures[W], set_velocity_futures)
                        , get_dependency<BACK>(nx_block, ny_block, nz_block, recv_futures[W], set_velocity_futures)
                        , get_dependency<BOTTOM_RIGHT>(nx_block, ny_block, nz_block, recv_futures[W], set_velocity_futures)
                        , get_dependency<BACK_BOTTOM>(nx_block, ny_block, nz_block, recv_futures[W], set_velocity_futures)

                    );

                 compute_fg_futures(nx_block, ny_block, nz_block) = calc_future;
            }

    send_boundaries_F(compute_fg_futures, step_);
    receive_boundaries_F(recv_futures, step_);

    send_boundaries_G(compute_fg_futures, step_);
    receive_boundaries_G(recv_futures, step_);

    send_boundaries_H(compute_fg_futures, step_);
    receive_boundaries_H(recv_futures, step_);

    for (std::size_t nz_block = 0; nz_block < c.num_z_blocks; ++nz_block)
        for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
            for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
            {
                hpx::shared_future<void> calc_future =
                    hpx::dataflow(
                        hpx::util::unwrapped(
                            hpx::util::bind(
                                &stencils<STENCIL_COMPUTE_RHS>::call,
                                boost::ref(rhs_data_),
                                boost::ref(data_[F]), boost::ref(data_[G]), boost::ref(data_[H]),
                                boost::ref(cell_type_data_),
                                boost::ref(fluid_cells_(nx_block, ny_block, nz_block)),
                                c.dx, c.dy, c.dz, dt
                            )
                        )
                        , compute_fg_futures(nx_block, ny_block, nz_block)
                        , get_dependency<LEFT>(nx_block, ny_block, nz_block, recv_futures[F], compute_fg_futures)
                        , get_dependency<FRONT>(nx_block, ny_block, nz_block, recv_futures[G], compute_fg_futures)
                        , get_dependency<BOTTOM>(nx_block, ny_block, nz_block, recv_futures[H], compute_fg_futures)
                    );

                 compute_rhs_futures(nx_block, ny_block, nz_block) = calc_future;
            }

    token.reset();
    for (std::size_t iter = 0; iter < c.iter_max; ++iter)
    {
         for (std::size_t nz_block = 0; nz_block < c.num_z_blocks; ++nz_block)
                for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
                    for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
                    {
                        hpx::shared_future<void> calc_future =
                            hpx::dataflow(
                                hpx::util::unwrapped(
                                    hpx::util::bind(
                                        &stencils<STENCIL_SET_P_OBSTACLE>::call,
                                        boost::ref(data_[P]),
                                        boost::ref(cell_type_data_),
                                        boost::ref(boundary_cells_(nx_block, ny_block, nz_block)),
                                        boost::ref(obstacle_cells_(nx_block, ny_block, nz_block)),
                                        token
                                    )
                                )
                                , compute_rhs_futures(nx_block, ny_block, nz_block)
                                , compute_res_futures(nx_block, ny_block, nz_block)
                            );

                         set_p_futures(nx_block, ny_block, nz_block) = calc_future;
                    }

            for (std::size_t nz_block = 0; nz_block < c.num_z_blocks; ++nz_block)
                for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
                    for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
                    {
                        hpx::shared_future<void> calc_future =
                            hpx::dataflow(
                                hpx::util::unwrapped(
                                    hpx::util::bind(
                                        &stencils<STENCIL_JACOBI>::call,
                                        boost::ref(data_[P]),
                                        boost::ref(rhs_data_),
                                        boost::ref(fluid_cells_(nx_block, ny_block, nz_block)),
                                        c.dx_sq, c.dy_sq, c.dz_sq, token
                                    )
                                )
                                , set_p_futures(nx_block, ny_block, nz_block)
                            );

                         solver_cycle_futures(nx_block, ny_block, nz_block) = calc_future;
                    }

            send_boundaries_P(solver_cycle_futures, step_ * c.iter_max + iter);
            receive_boundaries_P(recv_futures, step_ * c.iter_max + iter);

            for (std::size_t nz_block = 0; nz_block < c.num_z_blocks; ++nz_block)
                for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
                    for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
                    {
                        hpx::shared_future<Real> calc_future =
                            hpx::dataflow(
                                hpx::util::unwrapped(
                                    hpx::util::bind(
                                        &stencils<STENCIL_COMPUTE_RESIDUAL>::call,
                                        boost::ref(data_[P]),
                                        boost::ref(rhs_data_),
                                        boost::ref(fluid_cells_(nx_block, ny_block, nz_block)),
                                        c.dx_sq, c.dy_sq, c.dz_sq, token
                                    )
                                )
                                , solver_cycle_futures(nx_block, ny_block, nz_block)
                                , get_dependency<LEFT>(nx_block, ny_block, nz_block, recv_futures[P], solver_cycle_futures)
                                , get_dependency<RIGHT>(nx_block, ny_block, nz_block, recv_futures[P], solver_cycle_futures)
                                , get_dependency<BOTTOM>(nx_block, ny_block, nz_block, recv_futures[P], solver_cycle_futures)
                                , get_dependency<TOP>(nx_block, ny_block, nz_block, recv_futures[P], solver_cycle_futures)
                                , get_dependency<FRONT>(nx_block, ny_block, nz_block, recv_futures[P], solver_cycle_futures)
                                , get_dependency<BACK>(nx_block, ny_block, nz_block, recv_futures[P], solver_cycle_futures)
                            );

                         compute_res_futures(nx_block, ny_block, nz_block) = calc_future;
                    }

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

                if (c.rank == 0)
                {
                    hpx::future<std::vector<Real> > partial_residuals =
                        hpx::lcos::gather_here(residual_basename,
                                                std::move(local_residual),
                                                c.num_partitions, step_ * c.iter_max + iter, 0);

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
                            [dt, iter_int = iter, step = step_, t = t_, this](Real residual)
                            {
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

    local_max_uvs.clear();
    local_max_uvs.reserve(c.num_x_blocks * c.num_y_blocks * c.num_z_blocks);

    for (std::size_t nz_block = 0; nz_block < c.num_z_blocks; ++nz_block)
        for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
            for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
            {
                    local_max_uvs.push_back(hpx::dataflow(
                        hpx::util::unwrapped(
                            hpx::util::bind(
                                &stencils<STENCIL_UPDATE_VELOCITY>::call,
                                boost::ref(data_[U]), boost::ref(data_[V]), boost::ref(data_[W]),
                                boost::ref(data_[F]), boost::ref(data_[G]), boost::ref(data_[H]),
                                boost::ref(data_[P]),
                                boost::ref(cell_type_data_),
                                boost::ref(fluid_cells_(nx_block, ny_block, nz_block)),
                                dt, c.over_dx, c.over_dy, c.over_dz
                            )
                        )
                        , compute_res_futures(nx_block, ny_block, nz_block)
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

   // hpx::wait_all(compute_res_futures.data_);

    return local_max_uv.get();
}

}
}
}
