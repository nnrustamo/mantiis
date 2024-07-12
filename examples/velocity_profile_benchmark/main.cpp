#include <chrono>
using namespace std::chrono;
#include <omp.h>

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/ximgproc.hpp>

#include "src/lb2d.hpp"

int main()
{
    omp_set_num_threads(4);
    auto start = high_resolution_clock::now();
    // simulation conditions
    _GLOBAL_::Cl = 1.0e-9;
    _GLOBAL_::Fbody = 1.0e-7;
    _GLOBAL_::mfp = 6.2e-7;

    //  Shape
    int Nx = 64, Ny = 64;

    Shape shape(Nx, Ny);
    std::string folder = "input/kn = 10/";
    // shape.loadExistingModel(folder);
    shape.addHorizontalBoundary(0);
    shape.addHorizontalBoundary(Ny - 1);
    shape.calculateProperties(_GLOBAL_::Cl, _GLOBAL_::mfp);
    shape.writeToText(folder);

    // Grid
    lattice latt;
    Grid2D G(shape, latt, false);

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
    simulationDetails << "This was single grid run, the number of active cells: " << G.gridSize << std::endl;

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

