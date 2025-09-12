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
#include <thread>

#include "lb2d.hpp"

struct Settings 
{   
    bool load_existing_model = 1;
    bool dump_final_results = 1;
    int verbose = 1;
    int dump_every_timestep = 1;
    int dump_distributions_first_timestep = 0;
    double threshold = 1.0e-10;
    double mfp = _GLOBAL_::mfp;
    double Fbody = _GLOBAL_::Fbody;
    double dx = _GLOBAL_::Cl;
    double T = _GLOBAL_::T_phy;
    double P = _GLOBAL_::P_phy;
    
    int collision_method = 2;
    int wall_boundary = 2;
    std::string simulation_type = "transport";
};

template<typename T>
struct WriterThread 
{
    static void write(const std::string& fName, const std::vector<T>& data) 
    {
        IO::writeVectorToFile(fName, data);
    }
};

class MantiisApp
{

public:
    bool is_multiblock;
    int Nx;
    int Ny;
    std::string folder;
    int iters;

private:
    std::vector<std::thread> writer_threads; 

    bool load_existing_model = 1;
    bool dump_final_results = 1;
    double tol = 1.0e-10;
    int verbose = 1;
    int dump_every_timestep = 1;
    int dump_distributions_first_timestep = 1e6;

    int proc_id;
    int num_procs;
    int num_threads;

    int collision_method;
    int wall_boundary;
    std::string simulation_type = "transport";

    lattice latt;
    std::unique_ptr<Shape> shape;
    std::unique_ptr<Grid2D> grid;
    std::unique_ptr<LB2D> lb;

public:
    MantiisApp(int argc, char* argv[], const Settings& settings) 
    {   
        load_existing_model                = settings.load_existing_model;
        dump_final_results                 = settings.dump_final_results;
        verbose                            = settings.verbose;
        dump_every_timestep                = settings.dump_every_timestep;
        dump_distributions_first_timestep  = settings.dump_distributions_first_timestep;
        tol                                = settings.threshold;
 
        _GLOBAL_::mfp                      = settings.mfp;
        _GLOBAL_::Fbody                    = settings.Fbody;
        _GLOBAL_::Cl                       = settings.dx;
        _GLOBAL_::T_phy                    = settings.T;
        _GLOBAL_::P_phy                    = settings.P;
 
        collision_method                   = settings.collision_method;
        wall_boundary                      = settings.wall_boundary;
        simulation_type                    = settings.simulation_type;

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
        if (proc_id == 0 && load_existing_model) 
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
        switch (collision_method)
        {
            case 2:
            {
                lb->setCollision(&LB2D::MRTCollisionRegularized);
                break;
            }
            case 1:
            {
                lb->setCollision(&LB2D::MRTCollisionNonRegularized);
                break;
            }
            case 0:
            {
                lb->setCollision(&LB2D::SRTCollision);
                break;
            }
            default:
            {
                throw std::runtime_error(
                                "Incorrect keyword for collision: available options are\n"
                                "\t 2 for regularized MRT\n"
                                "\t 1 for simple MRT\n"
                                "\t 0 for simple BGK\n");

            }
        }

        switch (wall_boundary)
        {
            case 2:
            {
                lb->setWallBoundary(&LB2D::SRBBWall);
                break;
            }
            case 1:
            {
                lb->setWallBoundary(&LB2D::MDBBWall);
                break;
            }
            case 0:
            {
                lb->setWallBoundary(&LB2D::BBWall);
                break;
            }
            default:
            {
                throw std::runtime_error(
                                "Incorrect keyword for wall boundary method: available options are\n"
                                "\t 2 for SRBBWall\n"
                                "\t 1 for MDBBWall\n"
                                "\t 0 for BBWall\n");

            }
        }

        lb->setOpenBoundary(&LB2D::periodicBoundary);
    }

    void transport_simulation_loop()
    {
        lb->usq_old = lb->usq;
        lb->equilibrium();
        lb->evolutionStepOfMultiBlock(lb->grid.maxLevel);
        lb->calculateMacroscopicPorperties();
        lb->calculateVelocityDifference();
    }

    void adsorption_simulation_loop()
    {
        lb->rho_old = lb->rho;
        lb->calculateDensity();
        lb->getPseudoPotential();
        lb->calculateForces();
        lb->calculateVelocity();
        lb->equilibrium();
        lb->MRTCollisionNonRegularized();
        lb->applySourceTerm();
        lb->ftemp = lb->f; // save pre streaming
        lb->zeroDown();
        lb->stream();
        lb->BBWall();
        lb->periodicBoundary();
        lb->generalizedPeriodicBoundaryAdsorption();
        lb->calculateDensityDifference();
    }

    void coupled_simulation_loop()
    {
        lb->rho_old = lb->rho;
        lb->calculateDensity();
        lb->getPseudoPotential();
        lb->calculateForces();
        lb->calculateVelocity();
        lb->equilibrium();
        lb->MRTCollisionRegularized();
        lb->applySourceTerm();
        lb->ftemp = lb->f; // save pre streaming
        lb->zeroDown();
        lb->stream();
        lb->SRBBWall();
        lb->periodicBoundary();
        // lb->generalizedPeriodicBoundaryAdsorption();
        lb->calculateVelocityDifference();
    }

    void run()
    {   
        using namespace std::chrono;
        auto start = high_resolution_clock::now();
        std::string fName;

        auto run_loop = [&](auto simulation_step) 
        {
            for (int i = 0; i < iters; i++) 
            {
                simulation_step();

                if (proc_id == 0 && lb->t % verbose == 0)
                    std::cout << "Timestep: " << lb->t
                            << ", residual difference: " << lb->diff << "\n";

                if (lb->diff < tol)
                    break;

                if (lb->t < dump_distributions_first_timestep) 
                {
                    auto f_intermediate = lb->prepareDistributions();
                    fName = folder + "f_" + std::to_string(lb->t) + ".txt";
                    if (proc_id == 0)
                    {
                        // Join and remove finished threads to avoid resource leaks
                        for (auto it = writer_threads.begin(); it != writer_threads.end(); )
                        {
                            if (it->joinable())
                            {
                                it->join();
                                it = writer_threads.erase(it);
                            }
                            else
                            {
                                ++it;
                            }
                        }
                        // Launch a new thread for writing
                        writer_threads.emplace_back(&WriterThread<double>::write, fName, std::move(f_intermediate));
                    }
                }

                lb->t++;
            }
            // join all threads
            if (proc_id == 0)
            {
                for (auto& t : writer_threads)
                {
                    if (t.joinable())
                        t.join();
                }
                writer_threads.clear();
            }
        };

        if (simulation_type == "transport") 
        {
            run_loop([&] { transport_simulation_loop(); });
        }
        else if (simulation_type == "adsorption") 
        {
            run_loop([&] { adsorption_simulation_loop(); });
        }
        else if (simulation_type == "coupled") 
        {   
            std::cout<<"======================= Adsorption ONLY for initial condition =======================\n";
            auto fbody_temp = _GLOBAL_::Fbody;
            _GLOBAL_::Fbody = 0;
            auto iters_temp = iters;
            iters = 100;
            run_loop([&] 
                { adsorption_simulation_loop(); });
            iters = iters_temp;
            _GLOBAL_::Fbody = fbody_temp;
            std::fill(lb->ux.begin(), lb->ux.end(), 0);
            std::fill(lb->ux.begin(), lb->ux.end(), 0);
            std::fill(lb->rho.begin(), lb->rho.end(), _GLOBAL_::rhoin);
            std::fill(lb->feq.begin(), lb->feq.end(), 0);
             std::cout<<"======================= Coupling =======================\n";
            run_loop([&] { coupled_simulation_loop(); });
        }
        else 
        {
            throw std::runtime_error("Unknown simulation type: " + simulation_type);
        }

        std::vector<double> f_intermediate;
        if (dump_final_results)
        {   
            // lb->convertToPhysicalUnits();
            lb->completeSingleGridDomain();
            f_intermediate = lb->prepareDistributions();
        }

        if (proc_id == 0 && dump_final_results)
        {   
            IO::writeVectorToFile(folder + "f.txt", f_intermediate);

            IO::writeVectorToFile(folder + "ux.txt", lb->ux);
            IO::writeVectorToFile(folder + "uy.txt", lb->uy);
            IO::writeVectorToFile(folder + "rho.txt", lb->rho);
            IO::writeVectorToFile(folder + "convergence.txt", lb->diff_over_time);

            fName = folder + "SimulationDetails.txt";
            std::ofstream simulationDetails(fName);
            simulationDetails << "Domain dimension: " << shape->Nx << " x " << shape->Ny << "\n";
            simulationDetails << "Tolerance was set to " << tol << "\n";
            simulationDetails << "Number timesteps executed: " << lb->t << "\n";
            simulationDetails << "The number of active cells: " << grid->globalGridSize << "\n";
            auto end = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(end - start);
            std::cout << "Total run time (milliseconds) " << std::setprecision(15) << duration.count() << std::endl;
            simulationDetails << "Total run time (milliseconds): " << std::setprecision(15) << duration.count();
        }
    }
    ~MantiisApp(){MPI_Finalize();}
};
