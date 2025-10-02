#include <chrono>
#include "include/mantiisApp.hpp"

int main(int argc, char* argv[])
{   
    // sim settings
    Settings settings;

    settings.simulation_type = "coupled";
    settings.verbose = 10;
    settings.T = 333.0;
    settings.P = 10.0e6;
    settings.dx = 0.5e-11;

    auto little_focker = std::make_unique<MantiisApp>(argc, argv, settings);
    little_focker->init();
    little_focker->run();
    return 0;
}
