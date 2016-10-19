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
    cells_y_(c.cells_y_per_block * c.num_y_blocks + 2)
{

    if (c.verbose)
        std::cout << "Solver: blockwise Jacobi" << std::endl;

    if (c.verbose)
        std::cout << "Parellelization: custom grain size" << std::endl;

  //  if (c.verbose)
  //      std::cout << c << std::endl;

    step_ = 0;

    data_[U].resize(cells_x_, cells_y_, 0);
    data_[V].resize(cells_x_, cells_y_, 0);
    data_[F].resize(cells_x_, cells_y_, 0);
    data_[G].resize(cells_x_, cells_y_, 0);
    data_[P].resize(cells_x_, cells_y_, 0);
    data_[RHS].resize(cells_x_, cells_y_, 0);

    particles.reserve(c.num_fluid_cells * 16);

    Real dx_part = c.dx/4.;
    Real dy_part = c.dy/4.;

    cell_type_data_.resize(cells_x_,cells_y_);
    empty_marker_data_.resize(cells_x_, cells_y_);
    fluid_cells_.resize(c.num_x_blocks, c.num_y_blocks);
    obstacle_cells_.resize(c.num_x_blocks, c.num_y_blocks);

    for (std::size_t y = 1, ny_block = 0; y < cells_y_ - 1; y += c.cells_y_per_block, ++ny_block)
    {
        range_type y_range(y, std::min(y + c.cells_y_per_block, cells_y_ - 1));

        for (std::size_t x = 1, nx_block = 0; x < cells_x_ - 1; x += c.cells_x_per_block, ++nx_block)
        {
            range_type x_range(x, std::min(x + c.cells_x_per_block, cells_x_ - 1));

            for (std::size_t j = y_range.x; j < y_range.y; ++j)
                for (std::size_t i = x_range.x; i < x_range.y; ++i)
                {
                    cell_type_data_(i, j) = std::move(c.flag_grid[j * cells_x_ + i]);

                    auto empty_marker = c.empty_marker_grid[j * cells_x_ + i];

                    if (c.empty_marker_grid[j * cells_x_ + i - 1][is_empty_cell])
                        empty_marker[is_empty_cell_left] = 1;

                    if (c.empty_marker_grid[j * cells_x_ + i + 1][is_empty_cell])
                        empty_marker[is_empty_cell_right] = 1;

                    if (c.empty_marker_grid[(j - 1) * cells_x_ + i][is_empty_cell])
                        empty_marker[is_empty_cell_bottom] = 1;

                    if (c.empty_marker_grid[(j + 1) * cells_x_ + i][is_empty_cell])
                        empty_marker[is_empty_cell_top] = 1;

                    empty_marker_data_(i, j) = empty_marker;

                    if (cell_type_data_(i, j).test(is_fluid))
                    {
                        if (!empty_marker[is_empty_cell])
                        {
                            Real start_x = (i - 2) * c.dx + dx_part/2.;
                            Real start_y = (j - 2) * c.dy + dy_part/2.;

                            for (auto y = 0; y < 4; ++y)
                                for (auto x = 0; x < 4; ++x)
                                    particles.emplace_back(start_x + x * dx_part, start_y + y * dy_part);
                        }

                        fluid_cells_(nx_block, ny_block).emplace_back(i, j);
                    }
                    else if ((cell_type_data_(i, j).test(is_obstacle) && !cell_type_data_(i, j).test(is_boundary) && cell_type_data_(i, j).count() > 1)
                             || cell_type_data_(i, j).count() > 2)
                        obstacle_cells_(nx_block, ny_block).emplace_back(i, j);
                }

        }
    }

    for (auto y = cells_y_ - 1; y < cells_y_; --y)
    {
        for (auto x = 0; x < cells_x_; ++x)
            std::cout << empty_marker_data_(x, y).to_ulong() << " ";

        std::cout << std::endl;
    }
}

void partition_server::init()
{
    for (std::size_t var = 0; var < NUM_VARIABLES; ++var)
        data_[var].clear(0);

    step_ = 0;
    t_ = 0;
    next_out_ = -1e-12;
    outcount_ = 0;
    current = 1;
    last = 0;

    set_velocity_futures.resize(c.num_x_blocks, c.num_y_blocks);
    update_particle_futures.resize(c.num_x_blocks, c.num_y_blocks);
    compute_fg_futures.resize(c.num_x_blocks, c.num_y_blocks);
    compute_rhs_futures.resize(c.num_x_blocks, c.num_y_blocks);
    compute_res_futures.resize(c.num_x_blocks, c.num_y_blocks);

    for (auto& a : compute_res_futures)
        a = hpx::make_ready_future(0.);

    set_p_futures.resize(c.num_x_blocks, c.num_y_blocks);

    solver_cycle_futures[current].resize(c.num_x_blocks, c.num_y_blocks);
    solver_cycle_futures[last].resize(c.num_x_blocks, c.num_y_blocks);

    for (auto& a : solver_cycle_futures[last])
        a = hpx::make_ready_future();

    token.reset();
}

template<>
hpx::shared_future<void> partition_server::get_dependency<LEFT>(std::size_t idx_block, std::size_t idy_block, partition_data<hpx::shared_future<void> > calc_futures)
{
    if (idx_block == 0)
        return hpx::make_ready_future();

    return calc_futures(idx_block - 1, idy_block);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<RIGHT>(std::size_t idx_block, std::size_t idy_block, partition_data<hpx::shared_future<void> > calc_futures)
{
    if (idx_block == calc_futures.size_x_ - 1)
        return hpx::make_ready_future();

    return calc_futures(idx_block + 1, idy_block);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<BOTTOM>(std::size_t idx_block, std::size_t idy_block, partition_data<hpx::shared_future<void> > calc_futures)
{
    if (idy_block == 0)
        return hpx::make_ready_future();

    return calc_futures(idx_block, idy_block - 1);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<TOP>(std::size_t idx_block, std::size_t idy_block, partition_data<hpx::shared_future<void> > calc_futures)
{
    if (idy_block == calc_futures.size_y_ - 1)
        return hpx::make_ready_future();

    return calc_futures(idx_block, idy_block + 1);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<TOP_LEFT>(std::size_t idx_block, std::size_t idy_block, partition_data<hpx::shared_future<void> > calc_futures)
{
    if (idx_block == 0 || idy_block == calc_futures.size_y_ - 1)
        return hpx::make_ready_future();

    return calc_futures(idx_block - 1, idy_block + 1);
}

template<>
hpx::shared_future<void> partition_server::get_dependency<BOTTOM_RIGHT>(std::size_t idx_block, std::size_t idy_block, partition_data<hpx::shared_future<void> > calc_futures)
{
    if (idx_block == calc_futures.size_x_ - 1 || idy_block == 0)
        return hpx::make_ready_future();

    return calc_futures(idx_block + 1, idy_block - 1);
}

pair<Real> partition_server::do_timestep(Real dt)
{
    for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
    {
        for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
        {
           set_velocity_futures(nx_block, ny_block) =
                hpx::async(
                    hpx::util::bind(
                        &stencils<STENCIL_SET_VELOCITY>::call,
                        boost::ref(data_[U]), boost::ref(data_[V]),
                        boost::ref(cell_type_data_),
                        boost::ref(obstacle_cells_(nx_block, ny_block)),
                        c.bnd_condition
                    )
                );
        }
    }

    hpx::wait_all(set_velocity_futures.data_);

    for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
    {
        for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
        {
           update_particle_futures(nx_block, ny_block) =
                hpx::async(
                    hpx::util::bind(
                        &stencils<STENCIL_UPDATE_PARTICLE_POSITION>::call,
                        particles.begin(), particles.end(), boost::ref(data_[U]), boost::ref(data_[V]),
                        c.dx, c.dy, dt, c.x_length, c.y_length
                    )
                );
        }
    }

    if (c.vtk && next_out_ < t_)
    {
        next_out_ += c.delta_vec;

        if (c.verbose)
            std::cout << "Output in step " << step_ << " " << c.eps << " " << c.iter_max << token.was_cancelled() << std::endl;

        hpx::when_all(update_particle_futures.data_).then(
            hpx::launch::async,
            hpx::util::bind(
                &io::writer::write_vtk,
                boost::ref(data_[P]), boost::ref(data_[U]), boost::ref(data_[V]), boost::ref(cell_type_data_), boost::ref(particles),
                1, 1, c.i_max, c.j_max, c.dx, c.dx, outcount_++,
                0
            )
        ).wait();
    }

    for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
    {
        for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
        {
            compute_fg_futures(nx_block, ny_block) =
                hpx::dataflow(
                    hpx::util::unwrapped(
                        hpx::util::bind(
                            &stencils<STENCIL_COMPUTE_FG>::call,
                            boost::ref(data_[F]), boost::ref(data_[G]),
                            boost::ref(data_[U]), boost::ref(data_[V]),
                            boost::ref(cell_type_data_),
                            boost::ref(obstacle_cells_(nx_block, ny_block)),
                            boost::ref(fluid_cells_(nx_block, ny_block)),
                            c.re, c.gx, c.gy, c.dx, c.dy, dt, c.alpha
                        )
                    )
                    , set_velocity_futures(nx_block, ny_block)
                    , get_dependency<LEFT>(nx_block, ny_block, set_velocity_futures)
                    , get_dependency<LEFT>(nx_block, ny_block, set_velocity_futures)
                    , get_dependency<RIGHT>(nx_block, ny_block, set_velocity_futures)
                    , get_dependency<RIGHT>(nx_block, ny_block, set_velocity_futures)
                    , get_dependency<BOTTOM>(nx_block, ny_block, set_velocity_futures)
                    , get_dependency<BOTTOM>(nx_block, ny_block, set_velocity_futures)
                    , get_dependency<TOP>(nx_block, ny_block, set_velocity_futures)
                    , get_dependency<TOP>(nx_block, ny_block, set_velocity_futures)
                    , get_dependency<TOP_LEFT>(nx_block, ny_block, set_velocity_futures)
                    , get_dependency<BOTTOM_RIGHT>(nx_block, ny_block, set_velocity_futures)
                );
        }
    }

    for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
    {
        for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
        {
            compute_rhs_futures(nx_block, ny_block) =
                hpx::dataflow(
                    hpx::util::unwrapped(
                        hpx::util::bind(
                            &stencils<STENCIL_COMPUTE_RHS>::call,
                            boost::ref(data_[RHS]), boost::ref(data_[F]),
                            boost::ref(data_[G]), boost::ref(cell_type_data_),
                            boost::ref(fluid_cells_(nx_block, ny_block)),
                            c.dx, c.dy, dt
                        )
                    )
                    , compute_fg_futures(nx_block, ny_block)
                    , get_dependency<LEFT>(nx_block, ny_block, compute_fg_futures)
                    , get_dependency<BOTTOM>(nx_block, ny_block, compute_fg_futures)
                );
        }
    }

    token.reset();

    for (std::size_t iter = 0; iter < c.iter_max; ++iter)
    {
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
                                boost::ref(obstacle_cells_(nx_block, ny_block)),
                                token
                            )
                        )
                        , solver_cycle_futures[last](nx_block, ny_block)
                        , get_dependency<LEFT>(nx_block, ny_block, solver_cycle_futures[last])
                        , get_dependency<BOTTOM>(nx_block, ny_block, solver_cycle_futures[last])
                        , get_dependency<RIGHT>(nx_block, ny_block, solver_cycle_futures[last])
                        , get_dependency<TOP>(nx_block, ny_block, solver_cycle_futures[last])
                        , compute_res_futures(nx_block, ny_block)
                    );
            }
        }

        for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
        {
            for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
            {
                solver_cycle_futures[current](nx_block, ny_block) =
                    hpx::dataflow(
                        hpx::util::unwrapped(
                            hpx::util::bind(
                                &stencils<STENCIL_JACOBI>::call,
                                boost::ref(data_[P]), boost::ref(data_[RHS]),
                                boost::ref(fluid_cells_(nx_block, ny_block)),
                                c.dx_sq, c.dy_sq, token
                            )
                        )
                        , set_p_futures(nx_block, ny_block)
                        , solver_cycle_futures[last](nx_block, ny_block)
                        , get_dependency<LEFT>(nx_block, ny_block, solver_cycle_futures[last])
                        , get_dependency<BOTTOM>(nx_block, ny_block, solver_cycle_futures[last])
                        , get_dependency<RIGHT>(nx_block, ny_block, solver_cycle_futures[last])
                        , get_dependency<TOP>(nx_block, ny_block, solver_cycle_futures[last])
                        , compute_rhs_futures(nx_block, ny_block)
                    );
            }
        }

        for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
        {
            for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
            {
                compute_res_futures(nx_block, ny_block) =
                    hpx::dataflow(
                        hpx::util::unwrapped(
                            hpx::util::bind(
                                &stencils<STENCIL_COMPUTE_RESIDUAL>::call,
                                boost::ref(data_[P]), boost::ref(data_[RHS]),
                                boost::ref(fluid_cells_(nx_block, ny_block)), c.over_dx_sq, c.over_dy_sq,
                                token
                            )
                        )
                        , solver_cycle_futures[current](nx_block, ny_block)
                        , get_dependency<LEFT>(nx_block, ny_block,  solver_cycle_futures[current])
                        , get_dependency<RIGHT>(nx_block, ny_block, solver_cycle_futures[current])
                        , get_dependency<BOTTOM>(nx_block, ny_block, solver_cycle_futures[current])
                        , get_dependency<TOP>(nx_block, ny_block, solver_cycle_futures[current])
                    );
            }
        }

        std::swap(current, last);

        hpx::dataflow(
            hpx::util::unwrapped(
                [this, iter](std::vector<Real> residuals)
                {
                    if (!token.was_cancelled())
                    {
                        Real residual = 0;

                        for (std::size_t i = 0; i < residuals.size(); ++i)
                            residual += residuals[i];

                        residual = std::sqrt(residual)/c.num_fluid_cells;

                        if (residual < c.eps || iter == c.iter_max - 1)
                        {
                            std::cout << "Step: " << step_ << " Iterations: " << iter + 1 << " Residual: " << residual << std::endl;
                            token.cancel();
                        }
                    }

                }
            )
            , compute_res_futures.data_
        );
    }

    local_max_uvs.clear();
    local_max_uvs.reserve(c.num_x_blocks * c.num_y_blocks);

    for (std::size_t ny_block = 0; ny_block < c.num_y_blocks; ++ny_block)
    {
        for (std::size_t nx_block = 0; nx_block < c.num_x_blocks; ++nx_block)
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
                    , get_dependency<LEFT>(nx_block, ny_block, solver_cycle_futures[last])
                    , get_dependency<RIGHT>(nx_block, ny_block, solver_cycle_futures[last])
                    , get_dependency<BOTTOM>(nx_block, ny_block, solver_cycle_futures[last])
                    , get_dependency<TOP>(nx_block, ny_block, solver_cycle_futures[last])
                    , solver_cycle_futures[last](nx_block, ny_block)
                )
            );
        }
    }

    hpx::future<pair<Real> > local_max_uv =
            hpx::dataflow(
                hpx::util::unwrapped(
                    [](std::vector<pair<Real> > max_uvs)
                    -> pair<Real>
                    {
                        pair<Real>  max_uv(0);

                        for (std::size_t i = 0; i < max_uvs.size(); ++i)
                        {
                            max_uv.x = max_uvs[i].x > max_uv.x ? max_uvs[i].x : max_uv.x;
                            max_uv.y = max_uvs[i].y > max_uv.y ? max_uvs[i].y : max_uv.y;
                        }

                        return max_uv;
                    }
                )
                , local_max_uvs
            );

    t_ += dt;
    ++step_;

    return local_max_uv.get();
}

}
}
}
