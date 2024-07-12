#include <chrono>
using namespace std::chrono;
#include <omp.h>

#include "lb2d.hpp"

int main()
{
    omp_set_num_threads(14);
    auto start = high_resolution_clock::now();
    // simulation conditions
    _GLOBAL_::Cl = 1.0e-9;
    _GLOBAL_::Fbody = 1.0e-7;

    //  Shape
    int Nx = 512, Ny = 512;

    Shape shape(Nx, Ny);
    std::string folder = "in_out/";
    // shape.loadExistingModel(folder);
    
    shape.addHorizontalBoundary(0);
    shape.addHorizontalBoundary(Ny - 1);
    // shape.addRectangle(400, 600, 400, 600);
    shape.addCircle(20, 255, 255);
    shape.calculateProperties(_GLOBAL_::Cl, _GLOBAL_::mfp);
    shape.writeToText(folder);
    
    // Grid
    lattice latt;
    Grid2D G(shape, latt, true);

    // LB
    std::cout << "Initializing..." << std::endl;
    LB2D lb(Nx, Ny, latt, G, shape);
    lb.setCollision(&LB2D::MRTCollisionRegularized);
    lb.setOpenBoundary(&LB2D::periodicBoundary);
    lb.setWallBoundary(&LB2D::SRBBWall);

    double tol = 1.0e-6;
    int iter = 100000;
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
    std::cout << "Simulation run time (seconds) " << duration_sim.count() / 1000 << std::endl;

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);
    std::cout << "Total run time (seconds) " << duration.count() / 1000 << std::endl;

    // write domain and run details
    std::string fName = folder + "SimulationDetails.txt";
    std::ofstream simulationDetails(fName);
    simulationDetails << "Domain dimension: " << shape.Nx << " x " << shape.Ny << std::endl;
    simulationDetails << "Simulation run time in seconds: " << duration_sim.count() / 1000.0 << std::endl;
    simulationDetails << "Total run time in seconds: " << duration.count() / 1000.0 << std::endl;
    simulationDetails << "Tolerance was set to " << tol << std::endl;
    simulationDetails << "Number timesteps to converge: " << lb.t << std::endl;
    simulationDetails << "The number of active cells: " << G.gridSize << std::endl;

    // dump results
    lb.convertToPhysicalUnits();
    lb.ReconstructOriginalGrid();
    fName = folder + "ux.txt";
    IO::writeVectorToFile(fName, lb.ux);
    fName = folder + "uy.txt";
    IO::writeVectorToFile(fName, lb.uy);
    fName = folder + "rho.txt";
    IO::writeVectorToFile(fName, lb.rho);
    fName = folder + "convergence.txt";
    IO::writeVectorToFile(fName, lb.diff_over_time);

    return 0;
}

