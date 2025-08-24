#include <chrono>
#include "include/mantiisApp.hpp"

int main(int argc, char* argv[])
{   
    auto little_focker = std::make_unique<MantiisApp>(argc, argv);
    little_focker->init();
    little_focker->run();
    return 0;
}
