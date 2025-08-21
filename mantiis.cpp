#include <chrono>

#include "include/mantiisapp.hpp"

using namespace std::chrono;

int main(int argc, char* argv[])
{
    auto start = high_resolution_clock::now();
    
    MantiisApp mantiisapp(argc, argv);
    mantiisapp.init();
    mantiisapp.run();

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);
    std::cout << "Total run time (milliseconds) " << std::setprecision(15) << duration.count() << std::endl;

    return 0;
}
