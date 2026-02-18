#include <chrono>
#include "include/mantiisApp.hpp"

int main(int argc, char* argv[])
{   

    // sim settings
    Settings settings;

    settings.simulation_type = "transport";
    settings.verbose = 10;
    settings.T = 300.0;
    settings.P = 2.0e6;
    settings.dx = 5e-10;
    settings.Fbody = 1.0e-8;
    settings.dump_macroscopic_every_timestep = 1.0e6;
    settings.dump_distributions_first_timestep = -1;

    auto little_focker = std::make_unique<MantiisApp>(argc, argv, settings);
    little_focker->init();
    little_focker->run();
    return 0;
}
