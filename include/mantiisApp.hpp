/*
 * Copyright (c) January 2024
 *
 * Author: Nijat Rustamov
 * Email: nrustamo@uwyo.edu
 * Organization: University of Wyoming
 *
 * Academic Supervisor: Saman Aryana
 * Email: saryana@uwyo.edu
 * Organization: University of Wyoming
 *
 * This file is a part of Lattice Boltzmann Simulation Software
 * Proprietary Software - All Rights Reserved
 *
 * Unauthorized copying, modification, or distribution of this software,
 * or any portion of it is prohibited
 */

/*
    Main application class for Mantiis
*/

#pragma once

#include <omp.h>
#include <memory>
#include <chrono>

#include "lb2d.hpp"

class MantiisApp 
{
public:
    bool is_multiblock;
    int Nx;
    int Ny;
    std::string folder;

    int iters;
    double tol = 1.0e-10;
    int verbose = 1;
    int dump_every_timestep = 1;

    int proc_id;
    int num_procs;
    int num_threads;

    std::unique_ptr<MpiSession> mpi;
    lattice latt;
    std::unique_ptr<Shape> shape;
    std::unique_ptr<Grid2D> grid;
    std::unique_ptr<LB2D> lb;

public:
    MantiisApp(int argc, char* argv[]) 
    {
        if (argc < 6) 
        {
            std::cerr << "Usage: " << argv[0] << " <num_threads> <is_multiblock> <Nx> <input_folder/> <Nt>\n";
            std::exit(EXIT_FAILURE);
        }

        num_threads = std::stoi(argv[1]);
        is_multiblock = (std::string(argv[2]) == "true");
        Nx = std::stoi(argv[3]);
        Ny = Nx;
        folder = argv[4];
        iters = std::stoi(argv[5]);

        mantiis_parallel::launchMPI(proc_id, num_procs, num_threads);

        if (proc_id == 0) 
        {
            std::cout << "Mean free path: " << _GLOBAL_::mfp << "\n";
            std::cout << "Cl: " << _GLOBAL_::Cl << "\n";
            std::cout << "Fbody: " << _GLOBAL_::Fbody << "\n";
            std::cout << "Force conversion factor: " << std::setprecision(15) << _GLOBAL_::Cf << "\n";
            std::cout << "Input density: " << std::setprecision(15) << _GLOBAL_::rhoin << "\n";
        }
    }

    void init() 
    {   
        shape = std::make_unique<Shape>(Nx, Ny, _GLOBAL_::Cl);
        if (proc_id == 0) 
            shape->loadExistingModel(folder);  
        
        shape->MPI_broadcast(0);

        grid = std::make_unique<Grid2D>(*shape, latt, is_multiblock);
        if (proc_id == 0) 
        {
            grid->initialize(*shape);
            std::cout << "The number of active cells: " << grid->globalGridSize << "\n";
            std::cout << "Initializing Simulation...\n";
        }

        grid->MPI_broadcast(0);
        grid->MPI_distributeGridData();
        // grid->debugPrintVectors();

        lb = std::make_unique<LB2D>(Nx, Ny, latt, *grid, *shape);
        lb->initialize();
        lb->setCollision(&LB2D::MRTCollisionRegularized);
        lb->setOpenBoundary(&LB2D::periodicBoundary);
        lb->setWallBoundary(&LB2D::SRBBWall);
    }

    void run() 
    {   
        using namespace std::chrono;
        auto start = high_resolution_clock::now();

        std::string fName;
        for (int i = 0; i < iters; i++) 
        {
            lb->usq_old = lb->usq;
            lb->equilibrium();
            lb->evolutionStepOfMultiBlock(lb->grid.maxLevel);
            lb->calculateMacroscopicPorperties();
            lb->calculateVelocityDifference();

            if (verbose == 1 && proc_id == 0)
                std::cout << "Timestep: " << lb->t << ", error: " << lb->diff << "\n";

            if (lb->diff < tol)
                break;
            lb->t++;

            if (lb->t <= 100) 
            {
                // Optional: dump intermediate results
            }
        }

        // fName = folder + "SimulationDetails.txt";
        // std::ofstream simulationDetails(fName);
        // simulationDetails << "Domain dimension: " << shape->Nx << " x " << shape->Ny << "\n";
        // simulationDetails << "Tolerance was set to " << tol << "\n";
        // simulationDetails << "Number timesteps to converge: " << lb->t << "\n";
        // simulationDetails << "The number of active cells: " << grid->globalGridSize << "\n";

        // lb->ReconstructOriginalGrid();
        // IO::writeVectorToFile(folder + "ux.txt", lb->ux);
        // IO::writeVectorToFile(folder + "uy.txt", lb->uy);
        // IO::writeVectorToFile(folder + "rho.txt", lb->rho);
        // IO::writeVectorToFile(folder + "convergence.txt", lb->diff_over_time);
        
        if (proc_id == 0)
        {
            auto end = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(end - start);
            std::cout << "Total run time (milliseconds) " << std::setprecision(15) << duration.count() << std::endl;
        }
    }

    ~MantiisApp(){MPI_Finalize();}
};
