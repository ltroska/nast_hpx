#ifndef NAST_HPX_UTIL_CANCELLATION_TOKEN_HPP_
#define NAST_HPX_UTIL_CANCELLATION_TOKEN_HPP_

namespace nast_hpx { namespace util {
    
struct cancellation_token
{
private:
    typedef boost::atomic<bool> flag_type;
    std::shared_ptr<flag_type> was_cancelled_;

public:
    cancellation_token()
      : was_cancelled_(std::make_shared<flag_type>(false))
    {}

    bool was_cancelled() const HPX_NOEXCEPT
    {
        return was_cancelled_->load(boost::memory_order_relaxed);
    }

    void cancel() HPX_NOEXCEPT
    {
        was_cancelled_->store(true, boost::memory_order_relaxed);
    }
    
    void reset() HPX_NOEXCEPT
    {
        was_cancelled_->store(false, boost::memory_order_relaxed);
    }
};

}
}

#endif
