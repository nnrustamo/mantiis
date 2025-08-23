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

#include "lb2d.hpp"

class MantiisApp 
{
public:
    int num_threads;
    int num_procs;
    int proc_id;
    bool is_multiblock;
    int Nx;
    int Ny;
    std::string folder;

    int iters;
    double tol = 1.0e-10;
    int verbose = 1;
    int dump_every_timestep = 1;

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
        std::cout << "Process " << proc_id << " /" << num_procs << " launched\n";

        MPI_Barrier(MPI_COMM_WORLD);

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
        shape->loadExistingModel(folder);
        grid = std::make_unique<Grid2D>(*shape, latt, is_multiblock);
    
        if (proc_id == 0) 
        {
            std::cout << "The number of active cells: " << grid->globalGridSize << "\n";
            std::cout << "Initializing Simulation...\n";
        }

        grid->MPI_distributeGridData();
        // grid->debugPrintVectors();

        MPI_Barrier(MPI_COMM_WORLD);
        std::cout<<"Checkpoint 0\n";
        lb = std::make_unique<LB2D>(Nx, Ny, latt, *grid, *shape);
        std::cout<<"Checkpoint 1\n";
        std::cout<<"LB2D initialized with grid size: " << lb->grid.localGridSize << "\n";
        lb->setCollision(&LB2D::MRTCollisionRegularized);
        std::cout<<"Checkpoint 2\n";
        lb->setOpenBoundary(&LB2D::periodicBoundary);
        std::cout<<"Checkpoint 3\n";
        lb->setWallBoundary(&LB2D::SRBBWall);
        std::cout<<"Checkpoint 4\n";
        lb->initialize();
        std::cout<<"Checkpoint 5\n";
    }

    void run() 
    {
        std::string fName;
        for (int i = 0; i < iters; i++) 
        {
            lb->usq_old = lb->usq;
            lb->equilibrium();
            lb->evolutionStepOfMultiBlock(lb->grid.maxLevel);
            lb->calculateMacroscopicPorperties();
            lb->calculateVelocityDifference();

            if (verbose == 1)
                std::cout << "Timestep: " << lb->t << ", error: " << lb->diff << "\n";
            if (lb->diff < tol)
                break;
            lb->t++;

            if (lb->t <= 100) 
            {
                // Optional: dump intermediate results
            }
        }

        fName = folder + "SimulationDetails.txt";
        std::ofstream simulationDetails(fName);
        simulationDetails << "Domain dimension: " << shape->Nx << " x " << shape->Ny << "\n";
        simulationDetails << "Tolerance was set to " << tol << "\n";
        simulationDetails << "Number timesteps to converge: " << lb->t << "\n";
        simulationDetails << "The number of active cells: " << grid->globalGridSize << "\n";

        lb->ReconstructOriginalGrid();
        IO::writeVectorToFile(folder + "ux.txt", lb->ux);
        IO::writeVectorToFile(folder + "uy.txt", lb->uy);
        IO::writeVectorToFile(folder + "rho.txt", lb->rho);
        IO::writeVectorToFile(folder + "convergence.txt", lb->diff_over_time);
    }

    ~MantiisApp()
    {
        MPI_Finalize();
    }

};
