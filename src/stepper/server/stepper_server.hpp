#ifndef STEPPER_SERVER_STEPPER_HPP
#define STEPPER_SERVER_STEPPER_HPP

#include <hpx/include/components.hpp>
#include <hpx/lcos/local/receive_buffer.hpp>

#include "io/config.hpp"
#include "computation/parameters.hpp"
#include "computation/strategy.hpp"
#include "grid/partition.hpp"

#include <hpx/error.hpp>

namespace stepper { namespace server {

char const* stepper_basename = "/nast_hpx/stepper/";
char const* gather_basename = "/nast_hpx/gather/";

template <typename T, typename Mutex = hpx::lcos::local::spinlock>
    struct my_buffer
    {
    protected:
        typedef Mutex mutex_type;
        typedef hpx::lcos::local::promise<T> buffer_promise_type;

        struct entry_data
        {
        private:
            HPX_MOVABLE_BUT_NOT_COPYABLE(entry_data)

        public:
            entry_data()
              : can_be_deleted_(false)
            {}

            entry_data(entry_data && rhs)
              : promise_(std::move(rhs.promise_)),
                can_be_deleted_(rhs.can_be_deleted_)
            {}

            hpx::future<T> get_future()
            {
                return promise_.get_future();
            }

            template <typename Val>
            void set_value(Val && val)
            {
                promise_.set_value(std::forward<Val>(val));
            }

            buffer_promise_type promise_;
            bool can_be_deleted_;
        };

        typedef std::map<std::size_t, boost::shared_ptr<entry_data> >
            buffer_map_type;
        typedef typename buffer_map_type::iterator iterator;

        struct erase_on_exit
        {
            erase_on_exit(buffer_map_type& buffer_map, iterator it)
              : buffer_map_(buffer_map), it_(it)
            {}

            ~erase_on_exit()
            {
                buffer_map_.erase(it_);
            }

            buffer_map_type& buffer_map_;
            iterator it_;
        };

    private:
        HPX_MOVABLE_BUT_NOT_COPYABLE(my_buffer)

    public:
        my_buffer() {}

        my_buffer(my_buffer && other)
          : buffer_map_(std::move(other.buffer_map_))
        {}

        my_buffer& operator=(my_buffer && other)
        {
            if(this != &other)
            {
                buffer_map_ = std::move(other.buffer_map_);
            }
            return *this;
        }

        hpx::future<T> receive(std::size_t step)
        {
            boost::lock_guard<mutex_type> l(mtx_);

            iterator it = get_buffer_entry(step);
            HPX_ASSERT(it != buffer_map_.end());

            // if the value was already set we delete the entry after
            // retrieving the future
            if (it->second->can_be_deleted_)
            {
                erase_on_exit t(buffer_map_, it);
                return it->second->get_future();
            }

            // otherwise mark the entry as to be deleted once the value was set
            it->second->can_be_deleted_ = true;
            return it->second->get_future();
        }

        void store_received(std::size_t step, T && val)
        {
            boost::shared_ptr<entry_data> entry;

            {
                boost::lock_guard<mutex_type> l(mtx_);

                iterator it = get_buffer_entry(step);
                HPX_ASSERT(it != buffer_map_.end());

                entry = it->second;

                if (!entry->can_be_deleted_)
                {
                    // if the future was not retrieved yet mark the entry as
                    // to be deleted after it was be retrieved
                    entry->can_be_deleted_ = true;
                }
                else
                {
                    // if the future was already retrieved we can delete the
                    // entry now
                    buffer_map_.erase(it);
                }
            }

            // set value in promise, but only after the lock went out of scope
            entry->set_value(std::move(val));
        }

        void clear()
        {
            buffer_map_.clear();
        }

    protected:
        iterator get_buffer_entry(std::size_t step)
        {
            iterator it = buffer_map_.find(step);
            if (it == buffer_map_.end())
            {
                std::pair<iterator, bool> res =
                    buffer_map_.insert(
                        std::make_pair(step, boost::make_shared<entry_data>()));
                if (!res.second)
                {
                    HPX_THROW_EXCEPTION(hpx::error::invalid_status,
                        "base_receive_buffer::get_buffer_entry",
                        "couldn't insert a new entry into the receive buffer");
                }
                return res.first;
            }
            return it;
        }

    private:
        mutable mutex_type mtx_;
        buffer_map_type buffer_map_;
    };


struct HPX_COMPONENT_EXPORT stepper_server
    : hpx::components::component_base<stepper_server>
{
    public:
        stepper_server() {}
        stepper_server(uint num_localities);

        void setup(io::config cfg);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, setup, setup_action);

        void do_work();

        std::pair<RealType, RealType> do_timestep(uint step, RealType dt);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, do_timestep, do_timestep_action);

        void set_keep_running(uint iter, bool kr);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, set_keep_running, set_keep_running_action);

        void receive_p_action_(uint t, scalar_partition p, direction to_dir);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, receive_p_action_, receive_p_action);

        void receive_fg_action_(uint t, vector_partition fg, direction to_dir);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, receive_fg_action_, receive_fg_action);

        void receive_uv_action_(uint t, vector_partition uv, direction to_dir);
        HPX_DEFINE_COMPONENT_ACTION(stepper_server, receive_uv_action_, receive_uv_action);

    protected:
        template<typename T>
        void print_grid(std::vector<grid::partition<T> > const& grid, const std::string message = "") const;
        void write_vtk(uint step) const;

        uint get_index(uint k, uint l) const;

        void send_p_to_neighbor(uint t, scalar_partition p, direction dir);
        scalar_partition receive_p_from_neighbor(uint t, direction dir);
        void communicate_p_grid(uint iter);

        void send_fg_to_neighbor(uint t, vector_partition fg, direction dir);
        vector_partition receive_fg_from_neighbor(uint t, direction dir);
        void communicate_fg_grid(uint step);

        void send_uv_to_neighbor(uint t, vector_partition uv, direction dir);
        vector_partition receive_uv_from_neighbor(uint t, direction dir);
        void communicate_uv_grid(uint step);

    private:
        void initialize_parameters();
        void initialize_grids();
        void initialize_communication();

        RealType compute_new_dt(std::pair<RealType, RealType>) const;

        io::config c;
        computation::parameters params;
        computation::strategy* strategy;

        uint num_localities, num_localities_x, num_localities_y;
        std::vector<hpx::naming::id_type> localities;

        index_grid_type index_grid;
        vector_grid_type uv_grid, fg_grid;
        scalar_grid_type p_grid, rhs_grid;

        bool has_neighbor[NUM_DIRECTIONS];
        hpx::shared_future<hpx::id_type> neighbor_steppers_[NUM_DIRECTIONS];
        my_buffer<scalar_partition> p_recv_buffs_[NUM_DIRECTIONS];
        hpx::lcos::local::receive_buffer<vector_partition> fg_recv_buffs_[NUM_DIRECTIONS];
        hpx::lcos::local::receive_buffer<vector_partition> uv_recv_buffs_[NUM_DIRECTIONS];
        hpx::lcos::local::receive_buffer<bool> keep_running;

        scalar_partition scalar_dummy;
        vector_partition vector_dummy;
};

}//namespace server
}//namespace stepper

HPX_REGISTER_ACTION_DECLARATION(stepper::server::stepper_server::setup_action, stepper_server_setup_action);

#endif
