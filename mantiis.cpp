#include <chrono>
#include "include/mantiisApp.hpp"

int main(int argc, char* argv[])
{   
    // sim settings
    Settings settings;
    settings.dx = 1.0e-9;
    settings.Fbody = 0.00;
    settings.collision_method = 1;
    settings.wall_boundary = 0;
    settings.simulation_type = "adsorption";

    auto little_focker = std::make_unique<MantiisApp>(argc, argv, settings);
    little_focker->init();
    little_focker->run();
    return 0;
}
