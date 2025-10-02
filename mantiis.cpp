#include <chrono>
#include "include/mantiisApp.hpp"

int main(int argc, char* argv[])
{   
    double scaling_factor = 4.0;

    // sim settings
    Settings settings;

    settings.simulation_type = "transport";
    settings.verbose = 10;
    settings.T = 300.0;
    settings.P = 2.0e6;
    settings.dx = 1.0e-9 / scaling_factor;
    settings.Fbody = 1.0e-8 / scaling_factor;
    settings.dump_macroscopic_every_timestep = 1e6;

    auto little_focker = std::make_unique<MantiisApp>(argc, argv, settings);
    little_focker->init();
    little_focker->run();
    return 0;
}
