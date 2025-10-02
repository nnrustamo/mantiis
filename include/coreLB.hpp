/*
 * Copyright (c) January 2024
 *
 * Author: Nijat Rustamov
 * Organization: University of Wyoming
 * Email: nrustamo@uwyo.edu
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
    Core Lattice Boltzmann Method
*/

#pragma once

#include <vector>
#include <iostream>
#include <iomanip> // for std::setprecision
#include <cmath>
#include <map>
#include <string>

#include "linalg.h"
#include "stdafx.h"
#include "interpolation.h"
#include "ap.h"

#include "globalParameters.hpp"
#include "grid.hpp"
#include "comm.hpp"

// using namespace _GLOBAL_; // not the best way

class LB2D
{
private:
    const double F1 = 3.0; // in feq equation:  e/CS2
    const double F2 = 4.5; // in feq equation:  e/2*CS4
    const double F3 = 1.5; // in feq equation:  e/2*CS2
    const double CS = sqrt(1. / 3.);
    const double CS2 = 1. / 3.;

public:
    using MethodPtr = void (LB2D::*)(int);
    MethodPtr collisionMethod = nullptr;
    MethodPtr wallBoundaryMethod = nullptr;
    MethodPtr openBoundaryMethod = nullptr;

    // ================================================= Simulation domain
    int NX;
    int NY;
    double resolution = _GLOBAL_::Cl;
    int t = 0;

    // ================================================= D2Q9 lattice structure
    lattice latt;
    // Number of lattice directions
    const int NC = latt.e_x.size();
    // Weights of lattice directions
    const std::vector<double> W{4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};

    // ================================================= Large Data Vectors
    // Some of these following vectors may not be used
    // Depending on boundary types and simulation conditions
    std::vector<int> pore;
    std::vector<double> rho;
    std::vector<double> rho_old;
    std::vector<double> ux;
    std::vector<double> uy;
    std::vector<double> Fx;
    std::vector<double> Fy;
    std::vector<double> uysq;
    std::vector<double> uxsq;
    std::vector<double> usq;
    std::vector<double> usq_old;
    std::vector<double> usqMinusOld;

    std::vector<double> f;
    std::vector<double> fb; // pre-collision
    std::vector<double> ftemp;
    std::vector<double> feq;
    std::vector<double> fneq;
    std::vector<double> psi;
    std::vector<double> normal_x;
    std::vector<double> normal_y;
    std::vector<double> taus;

    // Index vectors
    std::vector<int64_t> streamID;
    std::vector<int64_t> streamLOC;
    std::vector<int> streamIC;
    std::vector<int64_t> boundaryID;
    std::vector<int> boundaryIC;

    std::vector<int64_t> periodicID;
    std::vector<int64_t> periodicI2;
    std::vector<int> periodicIC;

    std::vector<int64_t> generalizedPeriodicInlet;
    std::vector<int64_t> generalizedPeriodicInlet2;
    std::vector<int64_t> generalizedPeriodicOutlet;
    std::vector<int64_t> generalizedPeriodicOutlet2;

    std::vector<std::vector<int>> streamingForceIC;
    std::vector<std::vector<int64_t>> streamingForceI2;
    std::vector<std::vector<int>> bouncingForceIC;

    std::map<int64_t, int64_t> f_void_map;

    //
    double diff = 1.0;
    std::vector<double> diff_over_time;

    // grid
    Grid2D& grid;

    // mpi communcations
    std::unique_ptr<mantiis_parallel::MPICommunicationManager<double>> mpi_comm_manager;

public:
    LB2D(int &, int &, lattice &, Grid2D &, Shape &);

    void initialize();

    void equilibrium();

    void calculateMultigridTauS(Shape &);

    void calculateSingleGridTaus(Shape &);

    void getLocalMFP();

    void getPseudoPotential();

    void calculateForces();

    void applySourceTerm();

    void matMulLin(const std::vector<double> &, const std::vector<double> &, std::vector<double> &);

    void regularize();

    void SRTCollision(int lvl);

    void MRTCollisionNonRegularized(int);

    void MRTCollisionRegularized(int);

    void setCollision(MethodPtr);

    void setWallBoundary(MethodPtr);

    void setOpenBoundary(MethodPtr);

    void collide(int);

    void stream(int);

    void BBWall(int);

    void getBoundaryTypes();    // this function is deprecated

    void SRBBWall(int);

    void MDBBWall(int);

    void periodicBoundary(int);

    void generalizedPeriodicBoundaryAdsorption(int);

    void generalizedPeriodicBoundaryTransport(int);

    void applyBoundaryConditions(int);

    void calculateMacroscopicPorperties();

    void calculateDensity();

    void calculateVelocity();

    void interpolateBlockInterface(int, int);

    void zeroDown(int);

    void calculateVelocityDifference();

    void calculateDensityDifference();

    void simulate(double, int, int);

    void evolutionStepOfMultiBlock(int);

    void convertToPhysicalUnits();

    // functions to prepare variables to save
    std::vector<double> prepareUx();
    std::vector<double> prepareUy();
    std::vector<double> prepareRho();
    std::vector<double> prepareDistributions();

    // bring back the solid points in simulation domain
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> completeSingleGridDomain();
    void reconstructOriginalGrid();
};

// ========================================================================================
// Constructor
LB2D::LB2D(int &x, int &y, lattice &latt, Grid2D &G, Shape &shape) : NX(x), NY(y), latt(latt), grid(G)
{   
    // MPI communication indices
    std::vector<int64_t> comm_ids;
    std::vector<int64_t> comm_ids_ng;
    std::vector<int64_t> comm_ids_ic;
    std::vector<int64_t> comm_rank;

    // set up MPI communication indices
    mantiis_parallel::getMPICommunicationIDS(grid.startID, grid.endID,  grid.globalGridSize,
        grid.gridConnect, grid.bndTypes,
        comm_ids, comm_ids_ic, comm_ids_ng, comm_rank);    
    
    // for debugging
    // mantiis_parallel::printMPICommunicationIDs(comm_ids, comm_ids_ic, 
    //     comm_ids_ng, comm_rank);
    
    std::vector<mantiis_parallel::Instruction> send_instrcutions;
    for (int64_t i = 0; i < comm_ids.size(); i++) 
    {
        mantiis_parallel::Instruction instruction;
        instruction.id = comm_ids[i];
        instruction.ic = comm_ids_ic[i];
        instruction.ng = comm_ids_ng[i];
        instruction.peer = comm_rank[i];
        
        send_instrcutions.push_back(instruction);
    }

    int rank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    mpi_comm_manager = std::make_unique<mantiis_parallel::MPICommunicationManager<double>>(rank, send_instrcutions);
    MPI_Barrier(MPI_COMM_WORLD);
    // mpi_comm_manager->printCommunicationMaps();

    if (grid.isMultigrid)
        calculateMultigridTauS(shape);
    else
        calculateSingleGridTaus(shape);
    // Extract indices
    streamingForceIC.resize(grid.localGridSize);
    streamingForceI2.resize(grid.localGridSize);
    bouncingForceIC.resize(grid.localGridSize);

    // temporary vectors
    std::vector<int> vtemp1;
    std::vector<int> vtemp2;
    std::vector<int64_t> vtemp3;

    // inlet outlet 
    int ny_c;
    if (NY % 2 == 1)
        ny_c = (NY + 1) / 2 - 1;
    else
        ny_c = NY / 2 - 1;

    int i_in = ny_c;
    int j_in = 0;
    int i_out = ny_c;
    int j_out = NX - 1;

    for (int64_t i = 0; i < grid.localGridSize; i++)
    {
        for (int ic = 0; ic < NC; ic++)
        {   
            // check process bounds
            if (grid.gridConnect[i][ic] >= grid.startID && grid.gridConnect[i][ic] < grid.endID)
            {   
                int64_t neighbor = grid.gridConnect[i][ic] - grid.startID;
                bool is_x_periodic = (abs(grid.gridIJ[i][0] - grid.gridIJ[neighbor][0]) > static_cast<int>(pow(2, grid.maxLevel - 1)));
                bool is_y_periodic = (abs(grid.gridIJ[i][1] - grid.gridIJ[neighbor][1]) > static_cast<int>(pow(2, grid.maxLevel - 1)));

                if (!is_x_periodic && !is_y_periodic)
                {    
                    bool isBoundary = (grid.bndTypes[i] != -1); // has at least 1 boundary connection
                    bool isStream  = !isBoundary || _GLOBAL_::ifStream[grid.bndTypes[i]][ic];
                    if (isStream) // (-1) represents no boundary connections
                    {
                        streamID.push_back(i);
                        streamIC.push_back(ic);
                        streamLOC.push_back(neighbor);
                        vtemp1.push_back(ic);
                        vtemp3.push_back(neighbor);
                    }
                    else
                    {
                        boundaryID.push_back(i);
                        boundaryIC.push_back(ic);
                        vtemp2.push_back(ic);
                    }
                }
                else if (is_x_periodic)
                {
                    periodicID.push_back(i);
                    periodicI2.push_back(neighbor);
                    periodicIC.push_back(ic);
                    vtemp1.push_back(ic);
                    vtemp3.push_back(neighbor);
                }
            }
            else if (grid.gridConnect[i][ic] != -1 && (grid.gridConnect[i][ic] < grid.startID || grid.gridConnect[i][ic] >= grid.endID))
            {
                bool isBoundary = (grid.bndTypes[i] != -1); // has at least 1 boundary connection
                bool isStream = !isBoundary || _GLOBAL_::ifStream[grid.bndTypes[i]][ic];
                if (!isStream)
                {
                    boundaryID.push_back(i);
                    boundaryIC.push_back(ic);
                    vtemp2.push_back(ic);
                }
            }
            else if (grid.gridConnect[i][ic] == -1 && !_GLOBAL_::ifStream[grid.bndTypes[i]][ic])
            {
                boundaryID.push_back(i);
                boundaryIC.push_back(ic);
                vtemp2.push_back(ic);
            }
        }
        streamingForceIC[i] = vtemp1;
        bouncingForceIC[i] = vtemp2;
        streamingForceI2[i] = vtemp3;
        vtemp1.clear();
        vtemp2.clear();
        vtemp3.clear();

        // inlet outlet 
        if (grid.gridIJ[i][1] == i_in && grid.gridIJ[i][0] == j_in)
            generalizedPeriodicInlet.push_back(i);
        else if (grid.gridIJ[i][1] == i_out && grid.gridIJ[i][0] == j_out)
            generalizedPeriodicOutlet.push_back(i);
        else if (grid.gridIJ[i][1] == i_in && grid.gridIJ[i][0] == j_in + 1)
            generalizedPeriodicInlet2.push_back(i);
        else if (grid.gridIJ[i][1] == i_out && grid.gridIJ[i][0] == j_out - 1)
            generalizedPeriodicOutlet2.push_back(i);
    }

    if (rank == 0) 
    {
        int64_t counter = 0;
        for (int64_t i = 0; i < NX * NY; i++)
            if (!std::binary_search(grid.solidID.begin(), grid.solidID.end(), i)) 
            {
                f_void_map[i] = counter;
                counter++;
            }
    }
}

// Initialization
void LB2D::initialize()
{
    // Initialize distributions and other variables
    rho.resize(grid.localGridSize);
    rho_old.resize(grid.localGridSize);
    ux.resize(grid.localGridSize);
    uy.resize(grid.localGridSize);
    uxsq.resize(grid.localGridSize);
    uysq.resize(grid.localGridSize);
    usq.resize(grid.localGridSize);
    usq_old.resize(grid.localGridSize);
    usqMinusOld.resize(grid.localGridSize);
    f.resize(grid.localGridSize * NC);
    feq.resize(grid.localGridSize * NC);
    fb.resize(grid.localGridSize * NC);
    ftemp.resize(grid.localGridSize * NC);
    fneq.resize(grid.localGridSize * NC);
    psi.resize(grid.localGridSize);
    Fx.resize(grid.localGridSize);
    Fy.resize(grid.localGridSize);
    std::fill(rho.begin(), rho.end(), _GLOBAL_::rhoin);
    std::fill(ux.begin(), ux.end(), 0.0);
    std::fill(uy.begin(), uy.end(), 0.0);
    std::fill(f.begin(), f.end(), 0.0);

    int64_t k_ic;
#pragma omp parallel for default(shared) private(k_ic)
    for (int64_t k = 0; k < grid.localGridSize; k++)
    {
        k_ic = k * NC;
        for (int ic = 0; ic < NC; ic++)
        {
            f[k_ic + ic] = W[ic] * rho[k] * (1.0 + F1 * (latt.e_x[ic] * ux[k] + latt.e_y[ic] * uy[k]) + F2 * pow(latt.e_x[ic] * ux[k] + latt.e_y[ic] * uy[k], 2) - F3 * usq[k]);
        }
    }
    calculateMacroscopicPorperties();
}

// Equilibrium distribution
void LB2D::equilibrium()
{
    int64_t k_ic;
#pragma omp parallel for default(shared) private(k_ic)
    for (int64_t k = 0; k < grid.localGridSize; k++)
    {
        k_ic = k * NC;
        for (int ic = 0; ic < NC; ic++)
            feq[k_ic + ic] = W[ic] * rho[k] * (1.0 + F1 * (latt.e_x[ic] * ux[k] + latt.e_y[ic] * uy[k]) + F2 * pow(latt.e_x[ic] * ux[k] + latt.e_y[ic] * uy[k], 2) - F3 * usq[k]);
    }
}

void LB2D::calculateMultigridTauS(Shape &shape)
{
    // Calculate tau_s for each level
    taus.resize(grid.localGridSize);
    std::map<int, std::vector<std::vector<double>>> tauMap;
    for (int lvl = 0; lvl < grid.maxLevel; lvl++)
    {
        // Resize image
        double resizingFactor = pow(0.5, lvl);
        std::vector<std::vector<int>> resizedBinaryImage;
        utils::resizeBinaryImage(shape.domain, resizedBinaryImage, shape.Nx * resizingFactor, shape.Ny * resizingFactor);
        // Make a new shape
        Shape resizedShape(shape.Nx * resizingFactor, shape.Ny * resizingFactor, resolution / resizingFactor);
        resizedShape.addHorizontalBoundary(0);
        resizedShape.addHorizontalBoundary(resizedShape.Ny - 1);
        // Modify shape domain
        for (int i = 1; i < resizedShape.domain.size() - 1; i++)
        {
            for (int j = 0; j < resizedShape.domain[0].size(); j++)
            {
                resizedShape.domain[i][j] = resizedBinaryImage[i][j];
            }
        }
        // Calculate props
        resizedShape.calculateProperties(resolution / resizingFactor, _GLOBAL_::mfp);
        // Calculate taus
        std::vector<std::vector<double>> tau_level;
        tau_level.resize(resizedBinaryImage.size(), std::vector<double>(resizedBinaryImage[0].size(), 0.0));
        for (int i = 0; i < resizedShape.domain.size(); i++)
        {
            for (int j = 0; j < resizedShape.domain[0].size(); j++)
            {
                tau_level[i][j] = 0.5 + sqrt(6.0 / M_PI) * resizedShape.LocPore[i][j] * resizedShape.Kn[i][j] / (1 + 2 * resizedShape.Kn[i][j]);
            }
        }
        tauMap[lvl] = tau_level;
    }
    //
    int lvl, factor;
    std::vector<double> taus_temp(grid.globalGridSize);
    for (int64_t i = 0; i < grid.globalGridSize; i++)
    {
        lvl = grid.gridLevel[i];
        factor = pow(2, lvl - 1);
        std::vector<std::vector<double>> tau_level = tauMap[lvl - 1];
        taus_temp[i] = tau_level[grid.gridIJ[i][1] / factor][grid.gridIJ[i][0] / factor];
    }
    taus = mantiis_parallel::getSubDomainVector(taus_temp, grid.startID, grid.endID);
}

void LB2D::calculateSingleGridTaus(Shape &shape)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int64_t count = 0;
    taus.resize(grid.localGridSize);
    std::vector<double> taus_temp(grid.globalGridSize);
    for (int i = 0; i < shape.domain.size(); i++)
    {
        for (int j = 0; j < shape.domain[0].size(); j++)
        {
            if (shape.domain[i][j] == 1)
            {
                taus_temp[count] = 0.5 + sqrt(6.0 / M_PI) * shape.LocPore[i][j] * shape.Kn[i][j] / (1 + 2 * shape.Kn[i][j]);
                count++;
            }
        }
    }
    taus = mantiis_parallel::getSubDomainVector(taus_temp, grid.startID, grid.endID);
}

// local mean free path
void LB2D::getLocalMFP()
{
    for (int64_t i = 0; i < grid.localGridSize; i++)
    {
        _GLOBAL_::mfp = 1.0 / (sqrt(2.0) * _GLOBAL_::Crho_mol * rho[i] / _GLOBAL_::Mw_lu * _GLOBAL_::NA * M_PI * (2 * _GLOBAL_::d_mol / 2.0 + 2 * _GLOBAL_::d_mol / 2.0) * (2 * _GLOBAL_::d_mol / 2.0 + 2 * _GLOBAL_::d_mol / 2.0)); // pow(_GLOBAL_::d_mol, (double)2.0));
        grid.kn[i] = _GLOBAL_::mfp / (grid.locpore[i] * _GLOBAL_::Cl);
    }
}

// pseudo potential based on EOS
void LB2D::getPseudoPotential()
{
    double acc;
#pragma omp parallelgma omp parallel for default(shared) private(acc)
    for (int64_t i = 0; i < grid.localGridSize; i++)
    {
        double molar_density = rho[i] / _GLOBAL_::Mw_lu;
        acc = (molar_density * _GLOBAL_::R_lu * _GLOBAL_::T_lu / (1 - _GLOBAL_::b_lu * molar_density) - _GLOBAL_::a_alpha_lu * molar_density * molar_density / (1 + 2 * _GLOBAL_::b_lu * molar_density - _GLOBAL_::b_lu * _GLOBAL_::b_lu * molar_density * molar_density) - (1.0 / 3.0) * rho[i]) / (3.0 * _GLOBAL_::Gff);
        if (acc < 0)
        {
            std::cout << "pseudo potential is complex" << std::endl;
            psi[i] = acc;
        }
        else
            psi[i] = sqrt(acc);
    }
}

void LB2D::calculateForces()
{
    int64_t ic_curr, i_nex;
    double Fx_ff, Fy_ff, Fx_fw, Fy_fw;
    for (int64_t i = 0; i < grid.localGridSize; i++)
    {
        Fx_ff = 0;
        Fy_ff = 0;
        Fx_fw = 0;
        Fy_fw = 0;

        // the part where next node is not a wall
        for (int icx = 0; icx < streamingForceIC[i].size(); icx++)
        {
            ic_curr = streamingForceIC[i][icx];
            i_nex = streamingForceI2[i][icx];

            Fx_ff += (-6.0) * _GLOBAL_::Gff * psi[i] * _GLOBAL_::wf[ic_curr] * psi[i_nex] * latt.e_x[ic_curr];

            Fy_ff += (-6.0) * _GLOBAL_::Gff * psi[i] * _GLOBAL_::wf[ic_curr] * psi[i_nex] * latt.e_y[ic_curr];
        }

        // the part where next node is a wall
        for (int icx = 0; icx < bouncingForceIC[i].size(); icx++)
        {
            ic_curr = bouncingForceIC[i][icx];
            Fx_fw += (-1.0) * _GLOBAL_::Gfs * psi[i] * psi[i] * _GLOBAL_::wf[ic_curr] * latt.e_x[ic_curr];
            Fy_fw += (-1.0) * _GLOBAL_::Gfs * psi[i] * psi[i] * _GLOBAL_::wf[ic_curr] * latt.e_y[ic_curr];
        }

        Fx[i] = Fx_ff + Fx_fw + _GLOBAL_::Fbody;
        Fy[i] = Fy_ff + Fy_fw;
    }
}

void LB2D::applySourceTerm()
{
    int64_t k_ic;
    double cDotF;
    std::vector<double> uF(4);
    std::vector<double> ccs(4);
    std::vector<double> S(NC);
    std::vector<double> tempMat1(NC*NC);
    std::vector<double> tempMat2(NC);
    std::vector<double> tempMat3(NC);
    std::vector<double> tempMat4(NC);
    for (int64_t k = 0; k < grid.localGridSize; k++)
    {
        k_ic = k * NC;
        uF = {ux[k] * Fx[k], ux[k] * Fy[k],
              uy[k] * Fx[k], uy[k] * Fy[k]};

        for (int64_t ic = 0; ic < NC; ic++)
        {
            cDotF = latt.e_x[ic] * Fx[k] + latt.e_y[ic] * Fy[k];

            ccs = {(double)latt.e_x[ic] * latt.e_x[ic] - CS2, (double)latt.e_x[ic] * latt.e_y[ic],
                   (double)latt.e_y[ic] * latt.e_x[ic], (double)latt.e_y[ic] * latt.e_y[ic] - CS2};

            S[ic] = W[ic] * (cDotF / CS2 + (uF[0] * ccs[0] + uF[2] * ccs[1] + uF[1] * ccs[2] + uF[3] * ccs[3]) / (CS2 * CS2));
        }

        std::fill(tempMat1.begin(), tempMat1.end(), 0);
        std::fill(tempMat2.begin(), tempMat2.end(), 0);
        std::fill(tempMat3.begin(), tempMat3.end(), 0);
        std::fill(tempMat4.begin(), tempMat4.end(), 0);

        // relaxation parameters
        // double tau_s = 0.5 + sqrt((double)(6.0 / M_PI)) * (localPoreSize[k] - 1) * Kn[k] / (1.0 + 2.0 * Kn[k]);
        double tau_q = 0.5 + (3.0 + M_PI * (2.0 * taus[k] - 1.0) * (2.0 * taus[k] - 1.0) * _GLOBAL_::A2) / (8.0 * (2.0 * taus[k] - 1.0));
        double tau_d = 0.5 + (3.0 / 2.0) * sqrt((double)3.0) * 1.0 / pow((double)(1.0 / sqrt((double)(_GLOBAL_::mfp / _GLOBAL_::Cl * 1.0 / (1.0 + 2.0 * grid.kn[k]))) * 2.0), (double)2.0);

        _GLOBAL_::R = {
            1.0 / 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1.0 / 1.1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1.0 / 1.2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1.0 / tau_d, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1.0 / tau_q, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1.0 / tau_d, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1.0 / tau_q, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1.0 / taus[k], 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1.0 / taus[k]};

        // I - 1/2* R
        for (int u = 0; u < NC * NC; u++)
            tempMat1[u] = _GLOBAL_::I[u] - 0.5 * _GLOBAL_::R[u];

        // M x S
        matMulLin(_GLOBAL_::M, S, tempMat2);
        // (I - 1/2* R) x M x S
        matMulLin(tempMat1, tempMat2, tempMat3);
        // M-1 x (I - 1/2* R) x M x S
        matMulLin(_GLOBAL_::Minv, tempMat3, tempMat4);
        //
        for (int ic = 0; ic < NC; ic++)
            f[k_ic + ic] = f[k_ic + ic] + tempMat4[ic];
    }
}

// Collision step
// Function to multiply two linearized matrices
void LB2D::matMulLin(const std::vector<double> &mat1, const std::vector<double> &mat2, std::vector<double> &mat3)
{
    for (int u = 0; u < NC; u++)
    {
        for (int v = 0; v < NC; v++)
        {
            mat3[u] += mat1[NC * u + v] * mat2[v];
        }
    }
}

// Hermite regularization
void LB2D::regularize()
{
    std::vector<double> fneqtimescc(4, 0);
    std::vector<double> H2(4, 0);
    std::vector<double> ccr(4, 0);
    int64_t k_ic;
#pragma omp parallel for default(shared) private(k_ic) firstprivate(ccr, fneqtimescc, H2)
    for (int64_t k = 0; k < grid.localGridSize; k++)
    {
        std::fill(fneqtimescc.begin(), fneqtimescc.end(), 0);
        std::fill(H2.begin(), H2.end(), 0);
        std::fill(ccr.begin(), ccr.end(), 0);
        k_ic = k * NC;
        for (int ic = 0; ic < NC; ic++)
        {
            ccr = {(double)(latt.e_x[ic] * latt.e_x[ic]), (double)(latt.e_x[ic] * latt.e_y[ic]),
                   (double)(latt.e_x[ic] * latt.e_y[ic]), (double)(latt.e_y[ic] * latt.e_y[ic])};

            for (int u = 0; u < 4; u++)
                fneqtimescc[u] += (f[k_ic + ic] - feq[k_ic + ic]) * ccr[u];
        }

        for (int ic = 0; ic < NC; ic++)
        {
            H2 = {latt.e_x[ic] * latt.e_x[ic] / CS2 - 1, latt.e_x[ic] * latt.e_y[ic] / CS2, latt.e_x[ic] * latt.e_y[ic] / CS2, latt.e_y[ic] * latt.e_y[ic] / CS2 - 1};

            fneq[k_ic + ic] = W[ic] * (1. / (2. * CS2)) * (fneqtimescc[0] * H2[0] + fneqtimescc[1] * H2[1] + fneqtimescc[2] * H2[2] + fneqtimescc[3] * H2[3]);
        }
    }
}

// Collision functions
void LB2D::SRTCollision(int lvl = 1)
{
    int64_t k_ic;
#pragma omp parallel for default(shared) private(k_ic)
    for (int64_t i = 0; i < grid.localGridSize; i++)
    {
        if (grid.gridLevel[i] == lvl)
        {
            k_ic = i * NC;
            for (int ic = 0; ic < NC; ic++)
                f[k_ic + ic] = f[k_ic + ic] - 1.0 / _GLOBAL_::tau_SRT * (f[k_ic + ic] - feq[k_ic + ic]);
        }
    }
}

// MRT collision
void LB2D::MRTCollisionNonRegularized(int lvl = 1)
{
    std::vector<double> tempMat1(NC, 0);
    std::vector<double> tempMat2(NC, 0);
    std::vector<double> tempMat3(NC, 0);
    std::vector<double> tempMat4(NC, 0);
    int64_t k_ic;
#pragma omp parallel for default(shared) private(k_ic) firstprivate(_GLOBAL_::R, tempMat1, tempMat2, tempMat3, tempMat4)
    for (int64_t k = 0; k < grid.localGridSize; k++)
    {
        if (grid.gridLevel[k] == lvl)
        {
            k_ic = k * NC;
            // relaxation parameters
            // double tau_s = 0.5 + sqrt((double)(6.0 / M_PI)) * (localPoreSize[k] - 1) * Kn[k] / (1.0 + 2.0 * Kn[k]);
            double tau_q = 0.5 + (3.0 + M_PI * (2.0 * taus[k] - 1.0) * (2.0 * taus[k] - 1.0) * _GLOBAL_::A2) / (8.0 * (2.0 * taus[k] - 1.0));
            double tau_d = 0.5 + (3.0 / 2.0) * sqrt((double)3.0) * 1.0 / pow((double)(1.0 / sqrt((double)(_GLOBAL_::mfp / _GLOBAL_::Cl * 1.0 / (1.0 + 2.0 * grid.kn[k]))) * 2.0), (double)2.0);

            _GLOBAL_::R = {
                1.0 / 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1.0 / 1.1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1.0 / 1.2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1.0 / tau_d, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1.0 / tau_q, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1.0 / tau_d, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1.0 / tau_q, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1.0 / taus[k], 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1.0 / taus[k]};

            std::fill(tempMat1.begin(), tempMat1.end(), 0);
            std::fill(tempMat2.begin(), tempMat2.end(), 0);
            std::fill(tempMat3.begin(), tempMat3.end(), 0);
            std::fill(tempMat4.begin(), tempMat4.end(), 0);

            for (int ic = 0; ic < NC; ic++)
                tempMat1[ic] = f[k_ic + ic] - feq[k_ic + ic];

            // M x (f - feq)
            matMulLin(_GLOBAL_::M, tempMat1, tempMat2);
            // R x M x (f - feq)
            matMulLin(_GLOBAL_::R, tempMat2, tempMat3);
            // M-1 x R x M x (f - feq)
            matMulLin(_GLOBAL_::Minv, tempMat3, tempMat4);

            for (int ic = 0; ic < NC; ic++)
                f[k_ic + ic] = f[k_ic + ic] - tempMat4[ic];
        }
    }
}

void LB2D::MRTCollisionRegularized(int lvl = 1)
{
    regularize();
    std::vector<double> tempMat1(NC);
    std::vector<double> tempMat2(NC);
    std::vector<double> tempMat3(NC);
    std::vector<double> tempMat4(NC);
    int64_t k_ic;
#pragma omp parallel for default(shared) private(k_ic) firstprivate(_GLOBAL_::R, tempMat1, tempMat2, tempMat3, tempMat4)
    for (int64_t k = 0; k < grid.localGridSize; k++)
    {
        if (grid.gridLevel[k] == lvl)
        {
            k_ic = k * NC;
            // relaxation parameters
            double tau_q = 0.5 + (3.0 + M_PI * (2.0 * taus[k] - 1.0) * (2.0 * taus[k] - 1.0) * _GLOBAL_::A2) / (8.0 * (2.0 * taus[k] - 1.0));
            double tau_d = 1.0; //0.5 + (3.0 / 2.0) * sqrt((double)3.0) * 1.0 / pow((double)(1.0 / sqrt((double)(_GLOBAL_::mfp / _GLOBAL_::Cl * 1.0 / (1.0 + 2.0 * grid.kn[k]))) * 2.0), (double)2.0);

            _GLOBAL_::R = {
                1.0 / 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1.0 / 1.1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1.0 / 1.2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1.0 / tau_d, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1.0 / tau_q, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1.0 / tau_d, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1.0 / tau_q, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1.0 / taus[k], 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1.0 / taus[k]};

            std::fill(tempMat1.begin(), tempMat1.end(), 0);
            std::fill(tempMat2.begin(), tempMat2.end(), 0);
            std::fill(tempMat3.begin(), tempMat3.end(), 0);
            std::fill(tempMat4.begin(), tempMat4.end(), 0);

            for (int ic = 0; ic < NC; ic++)
                tempMat1[ic] = fneq[k_ic + ic];

            // M x (f - feq)
            matMulLin(_GLOBAL_::M, tempMat1, tempMat2);
            // R x M x (f - feq)
            matMulLin(_GLOBAL_::R, tempMat2, tempMat3);
            // M-1 x R x M x (f - feq)
            matMulLin(_GLOBAL_::Minv, tempMat3, tempMat4);

            for (int ic = 0; ic < NC; ic++)
                f[k_ic + ic] = feq[k_ic + ic] + fneq[k_ic + ic] - tempMat4[ic];
        }
    }
}

void LB2D::setCollision(MethodPtr method)
{
    collisionMethod = method;
}

void LB2D::setWallBoundary(MethodPtr method)
{
    wallBoundaryMethod = method;
}

void LB2D::setOpenBoundary(MethodPtr method)
{
    openBoundaryMethod = method;
}

void LB2D::collide(int lvl = 1)
{
    (this->*collisionMethod)(lvl);
}

void LB2D::applyBoundaryConditions(int lvl = 1)
{
    (this->*wallBoundaryMethod)(lvl);
    (this->*openBoundaryMethod)(lvl);

    // uncomment this to prepare data for scripts/debugger.py to analyze bugs
    // int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // std::ofstream fout("f_output_rank_" + std::to_string(rank) + ".txt");
    // for (int64_t i = 0; i < grid.localGridSize; i++)
    // {   
    //     for (int64_t ic = 0; ic < NC; ic++)
    //     {
    //         fout << "rank = " << rank << ", f[" << i + grid.startID << "][" << ic << "] = " << std::setprecision(16) << f[i * NC + ic] << std::endl;
    //     }
    // }
    // fout.close();
}

// Streaming steps
void LB2D::stream(int lvl = 1)
{
// std::fill(f.begin(), f.end(), 0);
#pragma omp parallel for default(shared)
    for (int64_t i = 0; i < streamID.size(); i++)
        if (grid.gridLevel[streamID[i]] == lvl)
            f[streamLOC[i] * NC + streamIC[i]] = ftemp[streamID[i] * NC + streamIC[i]];

    mpi_comm_manager->exchange(f, ftemp, MPI_COMM_WORLD);
    
    // uncomment this to prepare data for scripts/debugger.py to analyze bugs
    // int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // std::ofstream fout("f_output_rank_" + std::to_string(rank) + ".txt");
    // for (int64_t i = 0; i < grid.localGridSize; i++)
    // {   
    //     for (int64_t ic = 0; ic < NC; ic++)
    //     {
    //         fout << "rank = " << rank << ", f[" << i + grid.startID << "][" << ic << "] = " << std::setprecision(16) << f[i * NC + ic] << std::endl;
    //     }
    // }
    // fout.close();
}

// Simple halfway bounce back wall boundary
void LB2D::BBWall(int lvl = 1)
{
#pragma omp parallel for default(shared)
    for (int64_t i = 0; i < boundaryID.size(); i++)
        if (grid.gridLevel[boundaryID[i]] == lvl)
            f[boundaryID[i] * NC + latt.icOpp[boundaryIC[i]]] = ftemp[boundaryID[i] * NC + boundaryIC[i]];
}

void LB2D::getBoundaryTypes()
{
    /*
    This function is deprecated
    double cx_sum = 0, cy_sum = 0;
    bound.resize(grid.gridSize);
    std::vector<int64_t> connections;

    for (int64_t i = 0; i < grid.gridSize; i++)
    {
        cx_sum = 0;
        cy_sum = 0;
        connections = grid.gridConnect[i];
        for (int c = 0; c < connections.size(); c++)
        {
            if (connections[c] >= 0)
            {
                cx_sum += latt.e_x[c];
                cy_sum += latt.e_y[c];
            }
        }

        if ((cx_sum == 3 && cy_sum == 0) || (cx_sum == 2 && cy_sum == -1) || (cx_sum == 2 && cy_sum == 1) || (cx_sum == 1 && cy_sum == 0))
            bound[i] = 1;
        else if ((cx_sum == 0 && cy_sum == 3) || (cx_sum == 1 && cy_sum == 2) || (cx_sum == -1 && cy_sum == 2) || (cx_sum == 0 && cy_sum == 1))
            bound[i] = 2;
        else if ((cx_sum == -3 && cy_sum == 0) || (cx_sum == -2 && cy_sum == 1) || (cx_sum == -2 && cy_sum == -1) || (cx_sum == -1 && cy_sum == 0))
            bound[i] = 3;
        else if ((cx_sum == 0 && cy_sum == -3) || (cx_sum == 1 && cy_sum == -2) || (cx_sum == -1 && cy_sum == -2) || (cx_sum == 0 && cy_sum == -1))
            bound[i] = 4;
        else if (cx_sum == 1 && cy_sum == 1)
            bound[i] = 5;
        else if (cx_sum == -1 && cy_sum == 1)
            bound[i] = 6;
        else if (cx_sum == -1 && cy_sum == -1)
            bound[i] = 7;
        else if (cx_sum == 1 && cy_sum == -1)
            bound[i] = 8;
        else if ((cx_sum == 2 && cy_sum == 2) || (cx_sum == 1 && cy_sum == 3) || (cx_sum == 3 && cy_sum == 1))
            bound[i] = 9;
        else if ((cx_sum == -2 && cy_sum == 2) || (cx_sum == -1 && cy_sum == 3) || (cx_sum == -3 && cy_sum == 1))
            bound[i] = 10;
        else if ((cx_sum == -2 && cy_sum == -2) || (cx_sum == -3 && cy_sum == -1) || (cx_sum == -1 && cy_sum == -3))
            bound[i] = 11;
        else if ((cx_sum == 2 && cy_sum == -2) || (cx_sum == 1 && cy_sum == -3) || (cx_sum == 3 && cy_sum == -1))
            bound[i] = 12;
    }
    */
}

// Specular reflection + bounce back wall boundary
void LB2D::SRBBWall(int lvl = 1)
{
    int srOpp;
#pragma omp parallel for default(shared) private(srOpp, _GLOBAL_::r)
    for (int64_t i = 0; i < boundaryID.size(); i++)
        if (grid.gridLevel[boundaryID[i]] == lvl)
        {   
            _GLOBAL_::r = 0.4987*pow(grid.kn[boundaryID[i]], -0.131) - 0.25;
            // _GLOBAL_::r = 2.0 * A1 / (sqrt(6.0 / M_PI) + A1);
            // _GLOBAL_::r = 1/(1 + sqrt(M_PI/6.0) * A1 + (taus[boundaryID[i]] - 0.5) / (8 * pow(taus[boundaryID[i]] - 0.5, 2.0)));
            // std::cout<<"r value: "<<_GLOBAL_::r<<std::endl;
            srOpp = _GLOBAL_::icsr[grid.bndTypes[boundaryID[i]]][boundaryIC[i]];
            f[boundaryID[i] * NC + latt.icOpp[boundaryIC[i]]] += _GLOBAL_::r * ftemp[boundaryID[i] * NC + boundaryIC[i]];
            f[boundaryID[i] * NC + srOpp] += (1 - _GLOBAL_::r) * ftemp[boundaryID[i] * NC + boundaryIC[i]];
        }
}

// Maxwellian diffusion + bounce back wall boundary
void LB2D::MDBBWall(int lvl = 1)
{
    int64_t k1, k2;
    double Kno, Kdeno, normal_y, normal_x;
    for (int64_t i = 0; i < boundaryID.size(); i++)
    {   
        if (grid.gridLevel[boundaryID[i]] == lvl)
        {
            k1 = boundaryID[i];
            k2 = boundaryIC[i];
            Kdeno = 0.0;
            Kno = 0.0;
            normal_x = _GLOBAL_::xy_norm[grid.bndTypes[boundaryID[i]]][0];
            normal_y = _GLOBAL_::xy_norm[grid.bndTypes[boundaryID[i]]][1];
            for (int ic = 0; ic < NC; ic++)
            {
                if ((latt.e_y[ic] * normal_y + latt.e_x[ic] * normal_x) > 0)
                {
                    Kdeno += feq[k1 * NC + ic];
                }
                else if ((latt.e_y[ic] * normal_y + latt.e_x[ic] * normal_x) < 0)
                {
                    Kno += fb[k1 * NC + ic];
                }
            }
            _GLOBAL_::r = 0.3222*pow(grid.kn[boundaryID[i]], -0.217) - 0.25;
            // _GLOBAL_::r = 1/(1 + sqrt(M_PI/6.0) * A1 + (taus[boundaryID[i]] - 0.5) / (8 * pow(taus[boundaryID[i]] - 0.5, 2.0)));
            //std::cout<<"r value: "<<_GLOBAL_::r<<std::endl;
            f[k1 * NC + latt.icOpp[k2]] = _GLOBAL_::r * ftemp[k1 * NC + k2] + (1 - _GLOBAL_::r) * Kno / Kdeno * feq[k1 * NC + latt.icOpp[k2]];
        }
    }
}

void LB2D::periodicBoundary(int lvl = 1)
{
#pragma omp parallel for default(shared)
    for (int64_t i = 0; i < periodicID.size(); i++)
    {
        if (grid.gridLevel[periodicID[i]] == lvl)
            f[periodicI2[i] * NC + periodicIC[i]] = ftemp[periodicID[i] * NC + periodicIC[i]];
    }
}

// TBD
void LB2D::generalizedPeriodicBoundaryAdsorption(int lvl = 1)
{

#pragma omp parallel for default(shared)
    for (int64_t i = 0; i < generalizedPeriodicInlet.size(); i++)
    {
        int64_t k_1_in = generalizedPeriodicInlet[i];
        int64_t k_1_in_nc = k_1_in * NC;

        int64_t k_2_in = generalizedPeriodicInlet2[i];
        int64_t k_2_in_nc = k_2_in * NC;

        int64_t k_1_out = generalizedPeriodicOutlet[i];
        int64_t k_1_out_nc = k_1_out * NC;

        int64_t k_2_out = generalizedPeriodicOutlet2[i];
        int64_t k_2_out_nc = k_2_out * NC;


        // inlet
        f[k_1_in_nc + 1] = W[1] * _GLOBAL_::rhoin; // * (1.0 + F1 * ux[k_1_in] + F2 * uxsq[k_1_in] - F3 * usq[k_1_in]) + f[k_2_out_nc + 1] - feq[k_2_out_nc + 1];
        f[k_1_in_nc + 5] = W[2] * _GLOBAL_::rhoin; // * (1.0 + F1 * (ux[k_1_in] + uy[k_1_in]) + F2 * (ux[k_1_in] + uy[k_1_in]) * (ux[k_1_in] + uy[k_1_in]) - F3 * usq[k_1_in]) + f[k_2_out_nc + 5] - feq[k_2_out_nc + 5];
        f[k_1_in_nc + 8] = W[2] * _GLOBAL_::rhoin; // * (1.0 + F1 * (ux[k_1_in] - uy[k_1_in]) + F2 * (ux[k_1_in] - uy[k_1_in]) * (ux[k_1_in] - uy[k_1_in]) - F3 * usq[k_1_in]) + f[k_2_out_nc + 8] - feq[k_2_out_nc + 8];

        // outlet
        f[k_1_out_nc + 3] = W[1] * _GLOBAL_::rhoin; // * (1.0 - F1 * ux[k_1_out] + F2 * uxsq[k_1_out] - F3 * usq[k_1_out]) + f[k_2_in_nc + 3] - feq[k_2_in_nc + 3];
        f[k_1_out_nc + 6] = W[2] * _GLOBAL_::rhoin; // * (1.0 + F1 * (-ux[k_1_out] + uy[k_1_out]) + F2 * (-ux[k_1_out] + uy[k_1_out]) * (-ux[k_1_out] + uy[k_1_out]) - F3 * usq[k_1_out]) + f[k_2_in_nc + 6] - feq[k_2_in_nc + 6];
        f[k_1_out_nc + 7] = W[2] * _GLOBAL_::rhoin; // * (1.0 + F1 * (-ux[k_1_out] - uy[k_1_out]) + F2 * (-ux[k_1_out] - uy[k_1_out]) * (-ux[k_1_out] - uy[k_1_out]) - F3 * usq[k_1_out]) + f[k_2_in_nc + 7] - feq[k_2_in_nc + 7];
    }
}

// TBD
void LB2D::generalizedPeriodicBoundaryTransport(int lvl = 1)
{
}


//
void LB2D::calculateMacroscopicPorperties()
{
    double jx = 0, jy = 0;
    int64_t k_ic;
#pragma omp parallel for default(shared) reduction(+ : jx, jy) private(k_ic)
    for (int64_t i = 0; i < grid.localGridSize; i++)
    {
        jx = 0; jy = 0;

        rho[i] = 0;
        k_ic = i * NC;
        for (int64_t ic = 0; ic < NC; ic++)
        {

            rho[i] += f[k_ic + ic];

            jx += f[k_ic + ic] * latt.e_x[ic];

            jy += f[k_ic + ic] * latt.e_y[ic];
        }

        ux[i] = jx / rho[i] + _GLOBAL_::Fbody / (2.0 * rho[i]) * pow(2.0, grid.gridLevel[i] - 1);
        uy[i] = jy / rho[i];

        uxsq[i] = ux[i] * ux[i];
        uysq[i] = uy[i] * uy[i];
        usq[i] = uxsq[i] + uysq[i];

        // TODO
        if (ux[i] < 0)
        {
            // std::cout<<"negative velocity detected\n";
            // std::cout<<"Cell id: "<<grid.gridIJ[i][0]<<" "<<grid.gridIJ[i][1]<<std::endl;
            // exit(1);
        }
    }
}

void LB2D::calculateDensity()
{
    int64_t k_ic;
#pragma omp parallel for default(shared) private(k_ic)
    for (int64_t i = 0; i < grid.localGridSize; i++)
    {
        rho[i] = 0;
        k_ic = i * NC;
        for (int64_t ic = 0; ic < NC; ic++)
            rho[i] += f[k_ic + ic];
    }
}

void LB2D::calculateVelocity()
{
    double jx = 0, jy = 0;
    int64_t k_ic;
#pragma omp parallel for default(shared) reduction(+ : jx, jy) private(k_ic)
    for (int64_t i = 0; i < grid.localGridSize; i++)
    {
        jx = 0; jy = 0;
        k_ic = i * NC;
        for (int64_t ic = 0; ic < NC; ic++)
        {
            jx += f[k_ic + ic] * latt.e_x[ic];
            jy += f[k_ic + ic] * latt.e_y[ic];
        }

        ux[i] = jx / rho[i] + Fx[i] / (2.0 * rho[i]) * pow(2.0, grid.gridLevel[i] - 1);
        uy[i] = jy / rho[i] + Fy[i] / (2.0 * rho[i]) * pow(2.0, grid.gridLevel[i] - 1);;

        uxsq[i] = ux[i] * ux[i];
        uysq[i] = uy[i] * uy[i];
        usq[i] = uxsq[i] + uysq[i];

        // TODO
        if (ux[i] < 0)
        {
            // std::cout<<"negative velocity detected\n";
            // std::cout<<"Cell id: "<<grid.gridIJ[i][0]<<" "<<grid.gridIJ[i][1]<<std::endl;
            // exit(1);
        }
    }
}

void LB2D::interpolateBlockInterface(int l1, int l2)
{   
   //====================================================================== Coarse buffer cells
    int bufferType = l2 * 10 + l1;
    std::vector<int64_t> cells(5);
    double tau_coarse, tau_fine, avg_f;
#pragma omp parallel for default(shared) private(tau_fine, tau_coarse, avg_f, cells)
    for (int64_t l = 0; l < grid.coarseBuffer[bufferType].size(); l++)
    {
        cells = grid.coarseBuffer[bufferType][l];

        // calcualte coarse and fine relaxation
        tau_fine = 0;
        for (int ch_i = 1; ch_i < 5; ch_i++)
            tau_fine += taus[cells[ch_i]];
        tau_fine = tau_fine / 4;
        tau_coarse = taus[cells[0]];

        // average distributions
        for (int ic = 0; ic < NC; ic++)
        {
            avg_f = 0;
            for (int ch_i = 1; ch_i < 5; ch_i++)
                avg_f += f[cells[ch_i] * NC + ic];
            avg_f = avg_f / 4.0;

            // interpolation
            f[cells[0] * NC + ic] = feq[cells[0] * NC + ic] + 2 * tau_coarse / tau_fine * (avg_f - feq[cells[0] * NC + ic]);
        }
    }

    //====================================================================== Fine buffer cells
    // coordinate shift for centralization
    double parent_shift = pow(2.0, l2 - 2) - 0.5;
    double child_shift = (l1 == 1) ? 0 : pow(2.0, l1 - 2) - 0.5;
    bufferType = l1 * 10 + l2;
    std::vector<int64_t> parentConnections(9);

    alglib::real_1d_array x;
    alglib::real_1d_array y;
    alglib::real_1d_array Z;
    x.setlength(3);
    y.setlength(3);
    Z.setlength(9);

    std::vector<int> top_order = {7, 4, 8, 3, 0, 1, 6, 2, 5};
    std::vector<int> bottom_order = {6, 2, 5, 3, 0, 1, 7, 4, 8};
#pragma omp parallel for default(shared) private(tau_fine, tau_coarse) firstprivate(x, y, Z, cells, parentConnections)
    for (int64_t l = 0; l < grid.fineBuffer[bufferType].size(); l++)
    {
        cells = grid.fineBuffer[bufferType][l];
        parentConnections = grid.gridConnect[cells[0]];
        std::vector<double> X(9);
        std::vector<double> Y(9);
        // parent neighbor coordinates
        for (int ic = 0; ic < NC; ic++)
        {
            if (abs(grid.gridIJ[cells[0]][1] - grid.gridIJ[parentConnections[ic]][1]) > pow(2, grid.gridLevel[cells[0]] - 1) &&
                grid.gridIJ[cells[0]][1] - grid.gridIJ[parentConnections[ic]][1] < 0)
                X[ic] = grid.gridIJ[cells[0]][1] + parent_shift + 1 - 2;
            else if (abs(grid.gridIJ[cells[0]][1] - grid.gridIJ[parentConnections[ic]][1]) > pow(2, grid.gridLevel[cells[0]] - 1) &&
                        grid.gridIJ[cells[0]][1] - grid.gridIJ[parentConnections[ic]][1] > 0)
                X[ic] = grid.gridIJ[cells[0]][1] + parent_shift + 1 + 2;
            else
                X[ic] = grid.gridIJ[parentConnections[ic]][1] + parent_shift + 1;

            if (abs(grid.gridIJ[cells[0]][0] - grid.gridIJ[parentConnections[ic]][0]) > pow(2, grid.gridLevel[cells[0]] - 1) &&
                grid.gridIJ[cells[0]][0] - grid.gridIJ[parentConnections[ic]][0] < 0)
                Y[ic] = grid.gridIJ[cells[0]][0] + parent_shift + 1 - 2;

            else if (abs(grid.gridIJ[cells[0]][0] - grid.gridIJ[parentConnections[ic]][0]) > pow(2, grid.gridLevel[cells[0]] - 1) &&
                        grid.gridIJ[cells[0]][0] - grid.gridIJ[parentConnections[ic]][0] > 0)
                Y[ic] = grid.gridIJ[cells[0]][0] + parent_shift + 1 + 2;
            else
                Y[ic] = grid.gridIJ[parentConnections[ic]][0] + parent_shift + 1;
        }

        // sort and unique

        std::sort(X.begin(), X.end());
        auto last = std::unique(X.begin(), X.end());
        X.erase(last, X.end());
        for (int iX = 0; iX < X.size(); iX++)
        {
            // y[iX] = X[iX];

            if (grid.gridIJ[cells[0]][1] < NX / 2)
                y[iX] = X[iX];
            else
                y[iX] = X[X.size() - iX - 1];
        }

        std::sort(Y.begin(), Y.end());
        auto last1 = std::unique(Y.begin(), Y.end());
        Y.erase(last1, Y.end());
        for (int iX = 0; iX < Y.size(); iX++)
            x[iX] = Y[iX];
        X.clear();
        X.shrink_to_fit();
        Y.clear();
        Y.shrink_to_fit();
        for (int ic = 0; ic < NC; ic++)
        {
            // collect parent neighbor distributions
            for (int ic1 = 0; ic1 < NC; ic1++)
            {
                // Z[ic1] = f[parentConnections[top_order[ic1]] * NC + ic];

                if (grid.gridIJ[cells[0]][1] < NX / 2)
                    Z[ic1] = f[parentConnections[top_order[ic1]] * NC + ic];
                else
                    Z[ic1] = f[parentConnections[bottom_order[ic1]] * NC + ic];
            }
            // interpolate child cells
            // utils::BicubicInterpolator interp(X, Y, Z);
            // interp.solve();
            alglib::spline2dinterpolant s;
            alglib::spline2dbuildbicubicv(x, 3, y, 2, Z, 1, s);
            tau_coarse = taus[cells[0]];
            for (int ch_i = 1; ch_i < 5; ch_i++)
            {
                tau_fine = taus[cells[ch_i]];
                double f_interpolated = spline2dcalc(s, grid.gridIJ[cells[ch_i]][0] + child_shift + 1,
                                                        grid.gridIJ[cells[ch_i]][1] + child_shift + 1);
                f[cells[ch_i] * NC + ic] = feq[cells[ch_i] * NC + ic] +
                                            tau_fine / (2 * tau_coarse) *
                                                (f_interpolated - feq[cells[ch_i] * NC + ic]);
            }
        }
    }
}

void LB2D::zeroDown(int lvl = 1)
{
#pragma omp parallel for default(shared)
    for (int64_t i = 0; i < grid.localGridSize; i++)
    {
        if (grid.gridLevel[i] == lvl)
        {
            for (int ic = 0; ic < NC; ic++)
            {
                f[i * NC + ic] = 0;
            }
        }
    }
}

void LB2D::calculateVelocityDifference()
{   
    double sum1 = 0;
    double sum2 = 0;
#pragma omp parallel for default(shared) reduction(+ : sum1, sum2)
    for (int64_t i = 0; i < grid.localGridSize; i++)
    {
        usqMinusOld[i] = abs(sqrt(usq[i]) - sqrt(usq_old[i]));
        sum1 += usqMinusOld[i];
        sum2 += sqrt(usq[i]);
    }

    double sum1_cum, sum2_cum;

    MPI_Reduce(&sum1, &sum1_cum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sum2, &sum2_cum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    diff = sum1_cum / sum2_cum;
    MPI_Bcast(&diff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    diff_over_time.push_back(diff);
}

void LB2D::calculateDensityDifference()
{   
    double local_max = 0.0;

#pragma omp parallel for default(shared) reduction(max : local_max)
    for (int64_t i = 0; i < grid.localGridSize; i++)
    {
        auto rho_diff = std::abs(rho[i] - rho_old[i]) / rho_old[i];
        if (rho_diff > local_max)
            local_max = rho_diff;
    }
    if (local_max < 1.0e-15) local_max = 1.0e-5;
    MPI_Allreduce(&local_max, &diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    diff_over_time.push_back(diff);
}

// Main simulation loop
void LB2D::simulate(double tol = 1.0e-6, int iteration = 100000, int verbose = 0)
{
    initialize();
    t = 0;

    while (diff > tol && t < iteration)
    {
        usq_old = usq;
        equilibrium();
        collide();
        ftemp = f;
        zeroDown();
        stream(); // save post collision distribution function
        applyBoundaryConditions();
        calculateMacroscopicPorperties();
        calculateVelocityDifference();
        if (verbose == 1)
        {
            std::cout << "Timestep: " << t << ", error: " << diff << std::endl;
        }
        t++;
    }
}

void LB2D::evolutionStepOfMultiBlock(int lvl) // initial value of lvl (level) is MaxLevel
{
    // collide & stream as many times as the cells size
    int rep = pow(2, grid.maxLevel - lvl);
    for (int i = 0; i < rep; i++)
    {   
        fb = f; // save pre-collision
        collide(lvl);
        ftemp = f;  // save pre-streaming
        zeroDown(lvl);
        stream(lvl);
        applyBoundaryConditions(lvl);
    }

    // recursively call function itslef with lower levels
    if (lvl > 1)
        evolutionStepOfMultiBlock(lvl - 1);

    // interpolation start from the lowest level
    if (lvl != grid.maxLevel)
       interpolateBlockInterface(lvl, lvl + 1);
}

void LB2D::convertToPhysicalUnits()
{
    #pragma omp parallel for default(shared)
    for (int64_t i = 0; i < grid.localGridSize; i++)
    {
        ux[i] = ux[i] * _GLOBAL_::Cu;
        uy[i] = uy[i] * _GLOBAL_::Cu;
        rho[i] = rho[i] * _GLOBAL_::Crho;
    }
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> LB2D::completeSingleGridDomain()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Copy data
    std::vector<double> ux_copy(grid.localGridSize);
    std::vector<double> uy_copy(grid.localGridSize);
    std::vector<double> rho_copy(grid.localGridSize);

    #pragma omp parallel for default(shared)
    for (int64_t i = 0; i < grid.localGridSize; i++)
    {
        ux_copy[i] = ux[i];
        uy_copy[i] = uy[i];
        rho_copy[i] = rho[i];
    }

    // Delete original copy
    // ux.clear();
    // ux.shrink_to_fit();
    // uy.clear();
    // uy.shrink_to_fit();
    // rho.clear();
    // rho.shrink_to_fit();

    // Re-allocate memory
    // ux.resize(NX * NY);
    // std::fill(ux.begin(), ux.end(), 0.0);
    // uy.resize(NX * NY);
    // std::fill(uy.begin(), uy.end(), 0.0);
    // rho.resize(NX * NY);
    // std::fill(rho.begin(), rho.end(), 0.0);

    // Re-allocate memory
    std::vector<double> ux_complete(NX*NY, 0.0);
    std::vector<double> uy_complete(NX*NY, 0.0);
    std::vector<double> rho_complete(NX*NY, 0.0);

    // Gather all
    auto ux_gathered = mantiis_parallel::gatherToRoot(ux_copy);
    auto uy_gathered = mantiis_parallel::gatherToRoot(uy_copy);
    auto rho_gathered = mantiis_parallel::gatherToRoot(rho_copy);
    
    // Rewrite data with solid points
    if (rank == 0)
    {   
        std::sort(grid.solidID.begin(), grid.solidID.end());
        int64_t counter = 0;
        for (int64_t i = 0; i < NX * NY; i++)
        {
            if (std::binary_search(grid.solidID.begin(), grid.solidID.end(), i))
            {
                ux_complete[i] = 0.0;
                uy_complete[i] = 0.0;
                rho_complete[i] = 0.0;
            }
            else
            {
                ux_complete[i] = ux_gathered[counter];
                uy_complete[i] = uy_gathered[counter];
                rho_complete[i] = rho_gathered[counter];
                counter++;
            }
        }
    }
    return std::make_tuple(ux_complete, uy_complete, rho_complete);
}

void LB2D::reconstructOriginalGrid()
{
    // Copy data
    std::vector<double> ux_copy(grid.localGridSize);
    std::vector<double> uy_copy(grid.localGridSize);
    std::vector<double> rho_copy(grid.localGridSize);

    for (int64_t i = 0; i < grid.localGridSize; i++)
    {
        ux_copy[i] = ux[i];
        uy_copy[i] = uy[i];
        rho_copy[i] = rho[i];
    }

    // Delete original copy
    ux.clear();
    ux.shrink_to_fit();
    uy.clear();
    uy.shrink_to_fit();
    rho.clear();
    rho.shrink_to_fit();

    // Re-allocate memory
    ux.resize(NX * NY);
    std::fill(ux.begin(), ux.end(), 0);
    uy.resize(NX * NY);
    std::fill(uy.begin(), uy.end(), 0);
    rho.resize(NX * NY);
    std::fill(ux.begin(), ux.end(), 0);

    int i, j, dim, num_child_cells;

    double parent_shift;
    std::vector<int64_t> cellConnections(9);

    alglib::real_1d_array x;
    alglib::real_1d_array y;
    alglib::real_1d_array Z_ux;
    alglib::real_1d_array Z_uy;
    alglib::real_1d_array Z_rho;
    x.setlength(3);
    y.setlength(3);
    Z_ux.setlength(9);
    Z_uy.setlength(9);
    Z_rho.setlength(9);

    std::vector<int> top_order = {7, 4, 8, 3, 0, 1, 6, 2, 5};
    std::vector<int> bottom_order = {6, 2, 5, 3, 0, 1, 7, 4, 8};

    // Lowest level cells
    for (int64_t l = 0; l < grid.localGridSize; l++)
    {
        if (grid.gridLevel[l] == 1 && grid.gridIsBuffer[l] == 0)
        {
            i = grid.gridIJ[l][1];
            j = grid.gridIJ[l][0];
            ux[i * NX + j] = ux_copy[l];
            uy[i * NX + j] = uy_copy[l];
            rho[i * NX + j] = rho_copy[l];
        }
    }

    // interpolate higher level cells
    for (int64_t l = 0; l < grid.localGridSize; l++)
    {
        if (grid.gridLevel[l] > 1 && grid.gridIsBuffer[l] == 0)
        {
            // gather data for interpolation
            i = grid.gridIJ[l][1];
            j = grid.gridIJ[l][0];
            dim = grid.gridLevel[l];
            num_child_cells = pow(2, dim - 1);
            parent_shift = pow(2.0, static_cast<double>(dim) - 2) - 0.5;
            cellConnections = grid.gridConnect[l];

            std::vector<double> X(9);
            std::vector<double> Y(9);

            for (int ic = 0; ic < NC; ic++)
            {
                if (abs(grid.gridIJ[l][1] - grid.gridIJ[cellConnections[ic]][1]) > pow(2, grid.gridLevel[l] - 1) &&
                    grid.gridIJ[l][1] - grid.gridIJ[cellConnections[ic]][1] < 0)
                    X[ic] = grid.gridIJ[l][1] + parent_shift + 1 - 2;
                else if (abs(grid.gridIJ[l][1] - grid.gridIJ[cellConnections[ic]][1]) > pow(2, grid.gridLevel[l] - 1) &&
                         grid.gridIJ[l][1] - grid.gridIJ[cellConnections[ic]][1] > 0)
                    X[ic] = grid.gridIJ[l][1] + parent_shift + 1 + 2;
                else
                    X[ic] = grid.gridIJ[cellConnections[ic]][1] + parent_shift + 1;

                if (abs(grid.gridIJ[l][0] - grid.gridIJ[cellConnections[ic]][0]) > pow(2, grid.gridLevel[l] - 1) &&
                    grid.gridIJ[l][0] - grid.gridIJ[cellConnections[ic]][0] < 0)
                    Y[ic] = grid.gridIJ[l][0] + parent_shift + 1 - 2;

                else if (abs(grid.gridIJ[l][0] - grid.gridIJ[cellConnections[ic]][0]) > pow(2, grid.gridLevel[l] - 1) &&
                         grid.gridIJ[l][0] - grid.gridIJ[cellConnections[ic]][0] > 0)
                    Y[ic] = grid.gridIJ[l][0] + parent_shift + 1 + 2;
                else
                    Y[ic] = grid.gridIJ[cellConnections[ic]][0] + parent_shift + 1;

                if (grid.gridIJ[l][1] < NX / 2)
                {
                    Z_ux[ic] = ux_copy[cellConnections[top_order[ic]]];
                    Z_uy[ic] = uy_copy[cellConnections[top_order[ic]]];
                    Z_rho[ic] = rho_copy[cellConnections[top_order[ic]]];
                }
                else
                {
                    Z_ux[ic] = ux_copy[cellConnections[bottom_order[ic]]];
                    Z_uy[ic] = uy_copy[cellConnections[bottom_order[ic]]];
                    Z_rho[ic] = rho_copy[cellConnections[bottom_order[ic]]];
                }
            }

            // sort and unique
            std::sort(X.begin(), X.end());
            auto last = std::unique(X.begin(), X.end());
            X.erase(last, X.end());
            for (int iX = 0; iX < X.size(); iX++)
            {
                if (grid.gridIJ[l][1] < NX / 2)
                    y[iX] = X[iX];
                else
                    y[iX] = X[X.size() - iX - 1];
            }
            std::sort(Y.begin(), Y.end());
            auto last1 = std::unique(Y.begin(), Y.end());
            Y.erase(last1, Y.end());
            for (int iX = 0; iX < Y.size(); iX++)
                x[iX] = Y[iX];
            X.clear();
            X.shrink_to_fit();
            Y.clear();
            Y.shrink_to_fit();

            alglib::spline2dinterpolant s1;
            alglib::spline2dbuildbicubicv(x, 3, y, 2, Z_ux, 1, s1);
            alglib::spline2dinterpolant s2;
            alglib::spline2dbuildbicubicv(x, 3, y, 2, Z_uy, 1, s2);
            alglib::spline2dinterpolant s3;
            alglib::spline2dbuildbicubicv(x, 3, y, 2, Z_rho, 1, s3);

            // retrieve querry points
            std::vector<int> Xq(num_child_cells * num_child_cells);
            std::vector<int> Yq(num_child_cells * num_child_cells);
            std::vector<double> Zq_ux(num_child_cells * num_child_cells);
            std::vector<double> Zq_uy(num_child_cells * num_child_cells);
            std::vector<double> Zq_rho(num_child_cells * num_child_cells);

            // interpolate and assign
            double interpolated;
            for (int chi = i; chi < i + num_child_cells; chi++)
            {
                for (int chj = j; chj < j + num_child_cells; chj++)
                {
                    // interpolated = spline2dcalc(s1, chj + 1, chi + 1);
                    // ux[chi * NX + chj] = interpolated;
                    // interpolated = spline2dcalc(s2, chj + 1, chi + 1);
                    // uy[chi * NX + chj] = interpolated;
                    // interpolated = spline2dcalc(s3, chj + 1, chi + 1);
                    // rho[chi * NX + chj] = interpolated;

                    ux[chi * NX + chj] = ux_copy[l];
                    uy[chi * NX + chj] = uy_copy[l];
                    rho[chi * NX + chj] = rho_copy[l];
                }
            }
        }
    }
}

std::vector<double> LB2D::prepareUx()
{
    std::vector<double> ux_copy(NX*NY, 0.0);

    // fill in
    #pragma omp parallel for default(shared)
    for (int64_t i = 0; i < grid.gridID.size(); i++)
        ux_copy[grid.gridID[i]] = ux[i];

    return ux_copy;
}

std::vector<double> LB2D::prepareUy()
{
    std::vector<double> uy_copy(NX*NY, 0.0);

    // fill in
    #pragma omp parallel for default(shared)
    for (int64_t i = 0; i < grid.gridID.size(); i++)
        uy_copy[grid.gridID[i]] = uy[i];

    return uy_copy;
}

std::vector<double> LB2D::prepareRho()
{
    std::vector<double> rho_copy(NX*NY, 0.0);

    // fill in
    #pragma omp parallel for default(shared)
    for (int64_t i = 0; i < grid.gridID.size(); i++)
        rho_copy[grid.gridID[i]] = rho[i];

    return rho_copy;
}

std::vector<double> LB2D::prepareDistributions()
{   
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> f_copy(grid.localGridSize * NC, 0.0);
    // fill in
#pragma omp parallel for default(shared)
    for (int64_t i = 0; i < grid.gridID.size(); i++)
        for (int64_t ic = 0; ic < NC; ic ++)
            f_copy[i * NC + ic] = f[i * NC + ic];
    std::vector<double> f_gathered = mantiis_parallel::gatherToRoot(f_copy);
    std::vector<double> f_global(NX * NY * NC, 0.0);
    if (rank == 0) 
    {
        int64_t counter = 0;
        for (auto & pair: f_void_map)
        {
            for (int64_t ic = 0; ic < NC; ic++)
                f_global[pair.first * NC + ic] = f_gathered[pair.second * NC + ic];
        }
    }
    return f_global;
}
