#include <hpx/hpx_main.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/runtime/actions/plain_action.hpp>

#include "grid/partition_data.hpp"
#include "grid/partition.hpp"

void test(grid::partition const& left, grid::partition const& right) {
    hpx::cout << left.get_data(grid::left_partition).get();
}

HPX_PLAIN_ACTION(test, test_action);

int hpx_main()
{using namespace grid;


    std::vector<hpx::id_type> localities = hpx::find_all_localities();
    uint nl = localities.size();

    hpx::cout << "#localities = " << nl << "\n";

    partition par(hpx::find_here(), 3,3, 1);
    partition par2(hpx::find_here(), 3,3, 10);


    test_action act;
    hpx::async(act, hpx::find_here(), par, par2);

    // Initiate shutdown of the runtime system.
    return hpx::finalize();
}

int main()
{
    // Initialize and run HPX, this executes hpx_main above.
    return hpx::init();
}
