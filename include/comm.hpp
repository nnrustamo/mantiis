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

    std::tuple<int64_t, int64_t, int64_t, int64_t> calculateMPIGridPartition(int64_t gridSize)
    {
        int numProcs, rank;
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int64_t cellsPerProc = gridSize / numProcs;
        int64_t remainder = gridSize % numProcs;
        int64_t startIdx, endIdx, localGridSize;

        if (rank == numProcs - 1) 
        {
            startIdx = rank * cellsPerProc;
            localGridSize = cellsPerProc + remainder;
            endIdx = startIdx + localGridSize;
        } 
        else 
        {
            startIdx = rank * cellsPerProc;
            localGridSize = cellsPerProc;
            endIdx = startIdx + localGridSize;
        }

        return {startIdx, endIdx, localGridSize, cellsPerProc};
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

    void getMPICommunicationIDS(const int64_t startID, const int64_t endID, const int64_t cellsPerProc,
                                const std::vector<std::vector<int64_t>> &gridConnect,
                                std::vector<int64_t> &comm_ids, 
                                std::vector<int64_t> &comm_ids_ic,
                                std::vector<int64_t> &comm_ids_ng, 
                                std::vector<int64_t> &comm_rank)
    {
        
        int rank, numProcs;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

        for (int64_t i = 0; i < gridConnect.size(); i++)
            for (int ic = 0; ic < gridConnect[i].size(); ic++)
                if (gridConnect[i][ic] >= 0 && (gridConnect[i][ic] >= endID || gridConnect[i][ic] < startID))
                {
                    comm_ids.push_back(i);
                    comm_ids_ic.push_back(ic);
                    comm_ids_ng.push_back(gridConnect[i][ic] < startID ? cellsPerProc - (startID - gridConnect[i][ic]) : 
                                            gridConnect[i][ic] - endID);
                    comm_rank.push_back(gridConnect[i][ic] < startID ? (rank - 1 + numProcs) % numProcs : (rank + 1) % numProcs);
                }
    }
    
    void printMPICommunicationIDs(const std::vector<int64_t> &comm_ids,
                                  const std::vector<int64_t> &comm_ids_ic,
                                  const std::vector<int64_t> &comm_ids_ng,
                                  const std::vector<int64_t> &comm_rank)
    {
        int rank, numProcs;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
        
        for (int proc = 0; proc < numProcs; ++proc)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == proc)
            {
                std::cout << "Process Rank: " << rank << std::endl;
                std::cout << "Communication IDs:" << std::endl;
                for (size_t i = 0; i < comm_ids.size(); i++)
                {
                    std::cout << "  ID: " << comm_ids[i]
                              << ", IC: " << comm_ids_ic[i]
                              << ", NG: " << comm_ids_ng[i]
                              << ", Rank: " << comm_rank[i] << std::endl;
                }
                std::cout << std::flush;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    // Broadcast a 1D vector from root to all processes
    template <typename T>
    void broadcast_vector(std::vector<T> &vec, int root = 0)
    {   
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int64_t size = vec.size();

        MPI_Bcast(&size, 1, MPI_LONG_LONG, root, MPI_COMM_WORLD);

        if (rank != root)
            vec.resize(size);

        if constexpr (std::is_same_v<T, int>)
            MPI_Bcast(vec.data(), size, MPI_INT, root, MPI_COMM_WORLD);
        else if constexpr (std::is_same_v<T, int64_t>) 
            MPI_Bcast(vec.data(), size, MPI_LONG_LONG, root, MPI_COMM_WORLD);
        else if constexpr (std::is_same_v<T, double>)
            MPI_Bcast(vec.data(), size, MPI_DOUBLE, root, MPI_COMM_WORLD);
        else if constexpr (std::is_same_v<T, bool>)
        {
            std::vector<int> temp(size);
            if (rank == root)
                for (int64_t i = 0; i < size; i++) temp[i] = vec[i] ? 1 : 0;

            MPI_Bcast(temp.data(), size, MPI_INT, root, MPI_COMM_WORLD);

            if (rank != root)
                for (int64_t i = 0; i < size; i++) vec[i] = (temp[i] != 0);
        }
        else
            throw std::invalid_argument("Unsupported vector type for broadcast_vectors");
    }

    // Broadcast a 2D vector from root to all processes
    template <typename T>
    void broadcast_vector(std::vector<std::vector<T>> &vec2d, int root = 0)
    {   
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        int rows = vec2d.size();
        int cols = rows > 0 ? vec2d[0].size() : 0;
        MPI_Bcast(&rows, 1, MPI_INT, root, MPI_COMM_WORLD);
        MPI_Bcast(&cols, 1, MPI_INT, root, MPI_COMM_WORLD);

        if (rank != root)
        {
            vec2d.resize(rows);
            for (int i = 0; i < rows; ++i)
                vec2d[i].resize(cols);
        }
        
        std::vector<T> flat(rows * cols);
        if (root == 0)
        {
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    flat[i * cols + j] = vec2d[i][j];
        }

        if constexpr (std::is_same_v<T, int>)
            MPI_Bcast(flat.data(), rows * cols, MPI_INT, root, MPI_COMM_WORLD);
        else if constexpr (std::is_same_v<T, int64_t>)
            MPI_Bcast(flat.data(), rows * cols, MPI_LONG_LONG, root, MPI_COMM_WORLD);
        else if constexpr (std::is_same_v<T, double>)
            MPI_Bcast(flat.data(), rows * cols, MPI_DOUBLE, root, MPI_COMM_WORLD);
        else if constexpr (std::is_same_v<T, bool>)
        {
            std::vector<int> temp(rows * cols);
            if (root == 0)
                for (int i = 0; i < rows * cols; i++) temp[i] = flat[i];
            MPI_Bcast(temp.data(), rows * cols, MPI_INT, root, MPI_COMM_WORLD);
            if (root != 0)
                for (int i = 0; i < rows * cols; i++) flat[i] = temp[i];
        }
        else
            throw std::invalid_argument("Unsupported vector type for broadcast_vectors");
        
        if (rank != root)
        {
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    vec2d[i][j] = flat[i * cols + j];
        }
    }

} // namespace mantiis_parallel
