#include <chrono>
#include "include/mantiisApp.hpp"

int main(int argc, char* argv[])
{   
    // sim settings
    Settings settings;
    settings.dump_distributions_first_timestep = 1000;
    settings.dx = 1.0e-9;
    settings.simulation_type = "transport";
    settings.verbose = 10;

    auto little_focker = std::make_unique<MantiisApp>(argc, argv, settings);
    little_focker->init();
    little_focker->run();
    return 0;
}
