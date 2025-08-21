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
    Collection of MPI related functions
 */

#pragma once

#include <mpi.h>

namespace mantiis_parallel
{
    void launchMPI(int &proc_id, int &num_of_proc, int &num_thread)
    {
        omp_set_num_threads(num_thread);            
        MPI_Init(NULL, NULL);                       
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);    
        MPI_Comm_size(MPI_COMM_WORLD, &num_of_proc);
    }

    std::tuple<int64_t, int64_t, int64_t> calculateMPIGridPartition(int64_t gridSize)
    {
        int numProcs, rank;
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
        int64_t cellsPerProc = gridSize / numProcs;
        int64_t remainder = gridSize % numProcs;
    
        int64_t startIdx = rank * cellsPerProc + std::min(static_cast<int64_t>(rank), remainder);
        int64_t endIdx = startIdx + cellsPerProc + (rank < remainder ? 1 : 0);
        int64_t localGridSize = endIdx - startIdx;
    
        return {startIdx, endIdx, localGridSize};
    }

    template <typename T>
    T getSubDomainVector(const T &inputVector, int64_t startId, int64_t endId)
    {   
        if (startId >= inputVector.size() || endId > inputVector.size() || startId > endId)
            throw std::out_of_range("Invalid range for getSubDomainVector.");

        if constexpr (std::is_same_v<T, std::vector<typename T::value_type>>)
            return T(inputVector.begin() + startId, inputVector.begin() + endId);

        else if constexpr (std::is_same_v<T, std::vector<std::vector<typename T::value_type::value_type>>>)
            return T(inputVector.begin() + startId, inputVector.begin() + endId);
        else
            throw std::invalid_argument("Unsupported vector type.");
    }

} // namespace mantiis_parallel
