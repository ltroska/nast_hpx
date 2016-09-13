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

                            data_[U](i, j, k) = i + 100 * j + 10000 * k + 0.1 * hpx::get_locality_id();

                            cell_type_data_(i, j, k) = std::move(c.flag_grid[k * cells_x_ * cells_y_ + j * cells_x_ + i]);

                            if (cell_type_data_(i, j, k).test(is_fluid))
                                fluid_cells_(nx_block, ny_block, nz_block).emplace_back(i, j, k);
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
  //      data_[var].clear(0);

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
        recv_buffer_left_[P].valid_ = true;
        recv_buffer_left_[F].valid_ = true;
        recv_buffer_left_[G].valid_ = true;
        recv_buffer_left_[H].valid_ = true;
        recv_buffer_left_[U].valid_ = true;
        recv_buffer_left_[V].valid_ = true;
        recv_buffer_left_[W].valid_ = true;
    }

    if (!is_right_)
    {
        send_buffer_right_.dest_ = ids_[c.idz * c.num_partitions_x * c.num_partitions_y + c.idy * c.num_partitions_x + c.idx + 1];

        recv_buffer_right_[P].valid_ = true;
        recv_buffer_right_[U].valid_ = true;
        recv_buffer_right_[V].valid_ = true;
    }

    if (!is_bottom_)
    {
        send_buffer_bottom_.dest_ = ids_[(c.idz - 1) * c.num_partitions_x * c.num_partitions_y + c.idy * c.num_partitions_x + c.idx];

        recv_buffer_bottom_[U].valid_ = true;
    }

    if (!is_top_)
    {
        send_buffer_top_.dest_ = ids_[(c.idz + 1) * c.num_partitions_x * c.num_partitions_y + c.idy * c.num_partitions_x + c.idx];

        recv_buffer_top_[U].valid_ = true;
    }

    if (!is_front_)
    {
        send_buffer_front_.dest_ = ids_[c.idz * c.num_partitions_x * c.num_partitions_y + (c.idy - 1) * c.num_partitions_x + c.idx];

        recv_buffer_front_[U].valid_ = true;
    }

    if (!is_back_)
    {
        send_buffer_back_.dest_ = ids_[c.idz * c.num_partitions_x * c.num_partitions_y + (c.idy + 1) * c.num_partitions_x + c.idx];

        recv_buffer_back_[U].valid_ = true;
    }

    if (!is_back_ && !is_left_)
    {
        send_buffer_back_left_.dest_ = ids_[c.idz * c.num_partitions_x * c.num_partitions_y + (c.idy + 1) * c.num_partitions_x + c.idx - 1];

        recv_buffer_back_left_[U].valid_ = true;
    }

    if (!is_front_ && !is_right_)
    {
        send_buffer_front_right_.dest_ = ids_[c.idz * c.num_partitions_x * c.num_partitions_y + (c.idy - 1) * c.num_partitions_x + c.idx + 1];

        recv_buffer_front_right_[U].valid_ = true;
    }

    if (!is_bottom_ && !is_right_)
    {
        send_buffer_bottom_right_.dest_ = ids_[(c.idz - 1) * c.num_partitions_x * c.num_partitions_y + c.idy * c.num_partitions_x + c.idx + 1];

        recv_buffer_bottom_right_[U].valid_ = true;
    }

    if (!is_top_ && !is_left_)
    {
        send_buffer_top_left_.dest_ = ids_[(c.idz + 1) * c.num_partitions_x * c.num_partitions_y + c.idy * c.num_partitions_x + c.idx - 1];

        recv_buffer_top_left_[U].valid_ = true;
    }

    if (!is_back_ && !is_bottom_)
    {
        send_buffer_back_bottom_.dest_ = ids_[(c.idz - 1) * c.num_partitions_x * c.num_partitions_y + (c.idy + 1) * c.num_partitions_x + c.idx];

        recv_buffer_back_bottom_[U].valid_ = true;
    }

    if (!is_front_ && !is_top_)
    {
        send_buffer_front_top_.dest_ = ids_[(c.idz + 1) * c.num_partitions_x * c.num_partitions_y + (c.idy - 1) * c.num_partitions_x + c.idx];

        recv_buffer_front_top_[U].valid_ = true;
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
    recv_futures_P[current].resize(NUM_DIRECTIONS);
    recv_futures_P[last].resize(NUM_DIRECTIONS);

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


triple<Real> partition_server::do_timestep(Real dt)
{
   if (hpx::get_locality_id() == 0 || hpx::get_locality_id() == 1 || hpx::get_locality_id() == 2 || hpx::get_locality_id() == 4)
    {
        std::cout << "received " << std::endl;
        for (std::size_t z = data_[U].size_z_ - 1; z <= data_[U].size_z_ - 1; z--)
        {
            for (std::size_t y = data_[U].size_y_ - 1; y <= data_[U].size_y_ - 1; y--)
            {
                for (std::size_t x = 0; x < data_[U].size_x_; x++)
                    std::cout << data_[U](x, y, z) << " ";
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
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

    if (hpx::get_locality_id() == 1)
    {
        std::cout << " sending left" << std::endl;
        send_boundary<LEFT>(0, U, set_velocity_futures);
        std::cout << " done sending left " << std::endl;

        std::cout << " sending back left" << std::endl;
        send_boundary<BACK_LEFT>(0, U, set_velocity_futures);
        std::cout << " done sending back left " << std::endl;

        std::cout << " sending top left" << std::endl;
        send_boundary<TOP_LEFT>(0, U, set_velocity_futures);
        std::cout << " done sending top left " << std::endl;

            std::cout << " receiving left " << std::endl;

        receive_boundary<LEFT>(0, U, recv_futures);
        hpx::wait_all(recv_futures[U][LEFT].data_);
        std::cout << " done receiving left" << std::endl;

            std::cout << " receiving back left " << std::endl;

        receive_boundary<BACK_LEFT>(0, U, recv_futures);
        hpx::wait_all(recv_futures[U][BACK_LEFT].data_);
        std::cout << " done receiving back left" << std::endl;

            std::cout << " receiving top left " << std::endl;

        receive_boundary<TOP_LEFT>(0, U, recv_futures);
        hpx::wait_all(recv_futures[U][TOP_LEFT].data_);
        std::cout << " done receiving top left" << std::endl;
    }
    if (hpx::get_locality_id() == 0)
    {
                std::cout << " sending right" << std::endl;

        send_boundary<RIGHT>(0, U, set_velocity_futures);
                std::cout << " done sending right" << std::endl;

            std::cout << " sending top" << std::endl;

        send_boundary<TOP>(0, U, set_velocity_futures);
                std::cout << " done sending top" << std::endl;

        std::cout << " sending back" << std::endl;

        send_boundary<BACK>(0, U, set_velocity_futures);
                std::cout << " done sending back" << std::endl;

                        std::cout << " receiving right" << std::endl;

        receive_boundary<RIGHT>(0, U, recv_futures);
        hpx::wait_all(recv_futures[U][RIGHT].data_);
        std::cout << " done receiving right " << std::endl;

                   std::cout << " receiving TOP" << std::endl;

        receive_boundary<TOP>(0, U, recv_futures);
        hpx::wait_all(recv_futures[U][TOP].data_);
        std::cout << " done receiving TOP " << std::endl;

           std::cout << " receiving BACK" << std::endl;

        receive_boundary<BACK>(0, U, recv_futures);
        hpx::wait_all(recv_futures[U][BACK].data_);
        std::cout << " done receiving BACK " << std::endl;

    }

    if (hpx::get_locality_id() == 2)
    {
                std::cout << " sending front" << std::endl;

        send_boundary<FRONT>(0, U, set_velocity_futures);
                std::cout << " done sending front" << std::endl;

        std::cout << " sending front right" << std::endl;

        send_boundary<FRONT_RIGHT>(0, U, set_velocity_futures);
                std::cout << " done sending front right" << std::endl;

        std::cout << " sending front top" << std::endl;

        send_boundary<FRONT_TOP>(0, U, set_velocity_futures);
                std::cout << " done sending front front" << std::endl;

                        std::cout << " receiving front" << std::endl;

        receive_boundary<FRONT>(0, U, recv_futures);
        hpx::wait_all(recv_futures[U][FRONT].data_);
        std::cout << " done receiving front" << std::endl;

                        std::cout << " receiving front right" << std::endl;

        receive_boundary<FRONT_RIGHT>(0, U, recv_futures);
        hpx::wait_all(recv_futures[U][FRONT_RIGHT].data_);
        std::cout << " done receiving front right" << std::endl;

                        std::cout << " receiving front top" << std::endl;

        receive_boundary<FRONT_TOP>(0, U, recv_futures);
        hpx::wait_all(recv_futures[U][FRONT_TOP].data_);
        std::cout << " done receiving front top" << std::endl;

    }

     if (hpx::get_locality_id() == 4)
    {
                std::cout << " sending bottom" << std::endl;

        send_boundary<BOTTOM>(0, U, set_velocity_futures);
                std::cout << " done sending bottom" << std::endl;

            std::cout << " sending bottom right" << std::endl;

        send_boundary<BOTTOM_RIGHT>(0, U, set_velocity_futures);
                std::cout << " done sending bottom right" << std::endl;

            std::cout << " sending back bottom" << std::endl;

        send_boundary<BACK_BOTTOM>(0, U, set_velocity_futures);
                std::cout << " done sending back bottom" << std::endl;

                        std::cout << " receiving bottom " << std::endl;

        receive_boundary<BOTTOM>(0, U, recv_futures);
        hpx::wait_all(recv_futures[U][BOTTOM].data_);
        std::cout << " done receiving bottom" << std::endl;

                        std::cout << " receiving bottom right" << std::endl;

        receive_boundary<BOTTOM_RIGHT>(0, U, recv_futures);
        hpx::wait_all(recv_futures[U][BOTTOM_RIGHT].data_);
        std::cout << " done receiving bottom right" << std::endl;

                        std::cout << " receiving back bottom" << std::endl;

        receive_boundary<BACK_BOTTOM>(0, U, recv_futures);
        hpx::wait_all(recv_futures[U][BACK_BOTTOM].data_);
        std::cout << " done receiving back bottom" << std::endl;

    }

   if (hpx::get_locality_id() == 0 || hpx::get_locality_id() == 1 || hpx::get_locality_id() == 2 || hpx::get_locality_id() == 4)
    {
        std::cout << "received " << std::endl;
        for (std::size_t z = data_[U].size_z_ - 1; z <= data_[U].size_z_ - 1; z--)
        {
            for (std::size_t y = data_[U].size_y_ - 1; y <= data_[U].size_y_ - 1; y--)
            {
                for (std::size_t x = 0; x < data_[U].size_x_; x++)
                    std::cout << data_[U](x, y, z) << " ";
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }


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
