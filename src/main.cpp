#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

//#include "io/io.hpp"
//#include "stepper/stepper.hpp"
//#include "grid/types.hpp"

int hpx_main(int argc, char* argv[])
{
   /* if (argc != 2)
    {
        std::cerr << "Usage: ./main <input.xml>" << std::endl;
        exit(0);
    }*/

    auto c = hpx::dataflow(
        hpx::util::unwrapped(
            [](int a) -> int
            {return a;}),
                //hpx::make_ready_future(b)
            3
    );

   /* io::config config = io::read_config_from_file(argv[1]);

    stepper::stepper step;
    step.setup(config).wait();

  /*  std::cout << "create" << std::endl;
    auto a = vector_partition(hpx::find_here(), 3, 3);

    std::cout << "get data" << std::endl;
    auto data = a.get_data(CENTER).get();

    std::cout << "set data" << std::endl;
    data[0] = 4;

    std::cout << "cout data" << std::endl;
    std::cout << data << std::endl;

    std::cout << "get data and cout" << std::endl;
    std::cout << a.get_data(CENTER).get();

    std::cout << "new partition from data" << std::endl;
    a = vector_partition(hpx::find_here(), data);

    std::cout << "set data" << std::endl;
    data[0] = 3;

    std::cout << "cout data" << std::endl;
    std::cout << data << std::endl;

    std::cout << "get data and cout" << std::endl;
    std::cout << a.get_data(CENTER).get();

    std::cout << "cout data" << std::endl;
    std::cout << data << std::endl;

        std::cout << "done" << std::endl;*/

    return hpx::finalize();
}

int main(int argc, char* argv[])
{
    return hpx::init(argc, argv);
}
