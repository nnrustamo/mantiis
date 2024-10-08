#include <chrono>
#include <omp.h>
#include "src/lb2d.hpp"

using namespace std::chrono;

int main(int argc, char* argv[])
{
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <num_threads> <grid_type> <Nx> <input_folder/> <Nt> \n";
        std::cerr << "  <num_threads> : Number of OpenMP threads (e.g., 14)\n";
        std::cerr << "  <grid_type> : true or false (e.g., true)\n";
        std::cerr << "  <Nx> : Size of the grid in x-direction (e.g., 512)\n";
        std::cerr << "  <input_folder/> : Directory to read and write\n";
        std::cerr << "  <Nt> : Number of iterations\n";
        return 1;
    }

    // Parse command-line arguments
    int num_threads = std::stoi(argv[1]);
    bool initialize_grid = (std::string(argv[2]) == "true");

    omp_set_num_threads(num_threads);
    auto start = high_resolution_clock::now();
    
    // Simulation conditions
    // _GLOBAL_::T_phy = 300.0; // physical temperature, K;
    //_GLOBAL_::P_phy = 2.0e6; // physical pressure, Pa;
    _GLOBAL_::Cl = 1.0e-8;
    _GLOBAL_::Fbody = 1.0e-10;
    // _GLOBAL_::mfp  = 7.7308e-10;
    std::cout<<"Mean free path: "<<_GLOBAL_::mfp<<std::endl;
    std::cout<<"Cl: "<<_GLOBAL_::Cl<<std::endl;
    std::cout<<"Fbody: "<<_GLOBAL_::Fbody<<std::endl;
    std::cout<<"Focre conversion factor: "<<std::setprecision(15)<<_GLOBAL_::Cf<<std::endl;
    std::cout<<"Input density: "<<std::setprecision(15)<<_GLOBAL_::rhoin<<std::endl;

    // Shape
    int Nx = std::stoi(argv[3]);
    int Ny = Nx;
    Shape shape(Nx, Ny, _GLOBAL_::Cl);
    std::string folder = argv[4]; // "input_output/";
    shape.loadExistingModel(folder);
    // shape.addHorizontalBoundary(0);
    // shape.addHorizontalBoundary(Ny - 1);
    // shape.addRectangle(400, 600, 400, 600);
    // shape.addCircle(Nx/8, Nx/2, Ny/2);
    // shape.calculateProperties(_GLOBAL_::Cl, _GLOBAL_::mfp);
    // shape.writeToText(folder);

    // Grid
    lattice latt;
    Grid2D G(shape, latt, initialize_grid);

    // LB
    std::cout << "The number of active cells: " << G.gridSize << std::endl;
    std::cout << "Initializing Simulation..." << std::endl;
    LB2D lb(Nx, Ny, latt, G, shape);
    lb.setCollision(&LB2D::MRTCollisionRegularized);
    lb.setOpenBoundary(&LB2D::periodicBoundary);
    lb.setWallBoundary(&LB2D::SRBBWall);

    double tol = 1.0e-6;
    int iter = std::stoi(argv[5]);
    int verbose = 1;

    auto start_sim = high_resolution_clock::now();

    lb.initialize();
    for (int i = 0; i < iter; i++)
    {   
        lb.usq_old = lb.usq;
        lb.equilibrium();
        lb.evolutionStepOfMultiBlock(lb.grid.maxLevel);
        lb.calculateMacroscopicPorperties();
        lb.calculateVelocityDifference();
        if (verbose == 1)
        {
            std::cout << "Timestep: " << lb.t << ", error: " << lb.diff << std::endl;
        }
        if (lb.diff < tol)
        {
            break;
        }
        lb.t++;
    }

    auto end_sim = high_resolution_clock::now();
    auto duration_sim = duration_cast<milliseconds>(end_sim - start_sim);
    std::cout << "Simulation run time (milliseconds) " << std::setprecision(15) << duration_sim.count() << std::endl;

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);
    std::cout << "Total run time (milliseconds) " << std::setprecision(15) << duration.count() << std::endl;

    // Write domain and run details
    // std::string fName = folder + "SimulationDetails.txt";
    // std::ofstream simulationDetails(fName);
    // simulationDetails << "Domain dimension: " << shape.Nx << " x " << shape.Ny << std::endl;
    // simulationDetails << "Simulation run time in seconds: " << duration_sim.count() / 1000.0 << std::endl;
    // simulationDetails << "Total run time in seconds: " << duration.count() / 1000.0 << std::endl;
    // simulationDetails << "Tolerance was set to " << tol << std::endl;
    // simulationDetails << "Number timesteps to converge: " << lb.t << std::endl;
    // simulationDetails << "The number of active cells: " << G.gridSize << std::endl;

    // Dump results
    // lb.convertToPhysicalUnits();
    // lb.ReconstructOriginalGrid();
    // fName = folder + "ux.txt";
    // IO::writeVectorToFile(fName, lb.ux);
    // fName = folder + "uy.txt";
    // IO::writeVectorToFile(fName, lb.uy);
    // fName = folder + "rho.txt";
    // IO::writeVectorToFile(fName, lb.rho);
    // fName = folder + "convergence.txt";
    // IO::writeVectorToFile(fName, lb.diff_over_time);

    return 0;
}

