/*
 * Copyright (c) 2024 Nijat Rustamov
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
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/*
    Collection of MPI related functions
 */

#pragma once

#include <mpi.h>
#include <unistd.h>
#include "globalParameters.hpp"

constexpr int Q = 9; // hardcoded for D2Q9

namespace mantiis_parallel
{
    template<typename T> struct MpiType;
    template<> struct MpiType<int>       { static MPI_Datatype type() { return MPI_INT;       } };
    template<> struct MpiType<float>     { static MPI_Datatype type() { return MPI_FLOAT;     } };
    template<> struct MpiType<double>    { static MPI_Datatype type() { return MPI_DOUBLE;    } };
    template<> struct MpiType<int64_t>   { static MPI_Datatype type() { return MPI_LONG_LONG; } };

    struct Instruction 
    {
        int64_t id;   // row index
        int64_t ic;   // source col
        int64_t ng;   // dest col
        int64_t peer; // target rank
    };

    template<typename T>
    class MPICommunicationManager 
    {
    public:
        MPICommunicationManager(int64_t myrank,
                                const std::vector<Instruction>& send_instructions) : rank(myrank)
        {
            build(send_instructions);
        }

        void exchange(std::vector<T>& f,
              std::vector<T>& ftemp,
              MPI_Comm comm, int tagBase = 100) 
        {   
            int numProcs;
            MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
            MPI_Datatype dtype = MpiType<T>::type();
            int64_t count;

            std::vector<MPI_Request> requests;
            std::vector<std::vector<T>> all_send_buffers; // keep per-peer buffers alive

            // ========================= Post nonblocking sends =========================
            for (auto& peer_elem : send_map_id)
            {
                auto peer_rank = peer_elem.first;
                auto ids = peer_elem.second;
                auto ics = send_map_ic[peer_rank];

                std::vector<T> peer_send_buffer(ids.size()); // separate buffer per peer

                // pack data to send
                count = 0;
                for (size_t idx = 0; idx < ids.size(); ++idx)
                {
                    int64_t id = ids[idx];
                    int64_t ic = ics[idx];

                    if (count >= peer_send_buffer.size())
                        throw std::runtime_error("Send buffer overflow.");

                    peer_send_buffer[count] = ftemp[id * Q + ic];
                    count++;
                }

                MPI_Request req;
                MPI_Isend(peer_send_buffer.data(), peer_send_buffer.size(), dtype, peer_rank, tagBase, MPI_COMM_WORLD, &req);
                requests.push_back(req);

                // keep buffer alive until MPI_Waitall
                all_send_buffers.push_back(std::move(peer_send_buffer));
            }

            // ========================= Post nonblocking receives and unpack =========================
            for (auto& peer_elem : recv_map_id)
            {   
                auto peer_rank = peer_elem.first;
                recv_buffer.resize(peer_elem.second.size());

                MPI_Request req;
                MPI_Irecv(recv_buffer.data(), recv_buffer.size(), dtype, peer_rank, tagBase, MPI_COMM_WORLD, &req);
                MPI_Wait(&req, MPI_STATUS_IGNORE);

                // put sh*t together  
                auto ids = peer_elem.second;
                auto ics = recv_map_ic[peer_rank];
                count = 0;
                for (int64_t idx = 0; idx < ids.size(); ++idx)
                {
                    int64_t id = ids[idx];
                    int64_t ic = ics[idx];

                    if (count >= recv_buffer.size())
                        throw std::runtime_error("Receive buffer overflow.");
                    f[id * Q + ic] = recv_buffer[count];
                    count++;
                }
            }

            // ========================= Wait for all sends to complete =========================
            if (!requests.empty())
                MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
        }

        void printCommunicationMaps() const
        {
            MPI_Barrier(MPI_COMM_WORLD);
            sleep(2);
            std::cout << "==================== Send Map NG for "<<rank<<" with size "<<send_map_ng.size()<< " ====================" << std::endl;
            for (const auto& [key, values] : send_map_ng) 
            {
                std::cout << "Peer: " << key << " -> [ ";
                for (const auto& val : values) 
                {
                    std::cout << val << " ";
                }
                std::cout << "]" << std::endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
            sleep(2);
            std::cout << "==================== Send Map IC for "<<rank<<" with size "<<send_map_ic.size()<<" ====================" << std::endl;
            for (const auto& [key, values] : send_map_ic) 
            {
                std::cout << "Peer: " << key << " -> [ ";
                for (const auto& val : values) 
                {
                    std::cout << val << " ";
                }
                std::cout << "]" << std::endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
            sleep(2);
            std::cout << "==================== Recv Map ID for "<<rank<<" with size "<<recv_map_id.size()<<" ====================" << std::endl;
            for (const auto& [key, values] : recv_map_id) 
            {
                std::cout << "Peer: " << key << " -> [ ";
                for (const auto& val : values) 
                {
                    std::cout << val << " ";
                }
                std::cout << "]" << std::endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
            sleep(2);
            std::cout << "==================== Recv Map IC for "<<rank<<" with size "<<recv_map_ic.size()<<" ====================" << std::endl;
            for (const auto& [key, values] : recv_map_ic) 
            {
                std::cout << "Peer: " << key << " -> [ ";
                for (const auto& val : values) 
                {
                    std::cout << val << " ";
                }
                std::cout << "]" << std::endl;
            }

            MPI_Barrier(MPI_COMM_WORLD);
        }

    private:
        int64_t rank;
        std::unordered_map<int64_t, std::vector<int64_t>> send_map_id;
        std::unordered_map<int64_t, std::vector<int64_t>> send_map_ng;
        std::unordered_map<int64_t, std::vector<int64_t>> send_map_ic;
        std::unordered_map<int64_t, std::vector<int64_t>> recv_map_id;
        std::unordered_map<int64_t, std::vector<int64_t>> recv_map_ic;
        std::vector<T> send_buffer;
        std::vector<T> recv_buffer;
        
        void build(const std::vector<Instruction>& send_instructions)
        {      
            for (auto& in : send_instructions)
            {
                send_map_ng[in.peer].push_back(in.ng);
                send_map_ic[in.peer].push_back(in.ic);
                send_map_id[in.peer].push_back(in.id);
            }

            int tag = 0;

            // Persistent send buffers
            std::map<int, std::vector<int64_t>> send_buffers_ng;
            std::map<int, std::vector<int64_t>> send_buffers_ic;
            std::vector<MPI_Request> requests;
            
            for (auto& peer_send_data : send_map_ng)
            {
                int peer_rank = peer_send_data.first;

                // prepare send buffers
                send_buffers_ng[peer_rank] = peer_send_data.second;
                send_buffers_ic[peer_rank] = send_map_ic[peer_rank];

                MPI_Request req_ng, req_ic;

                // ========================= Send =========================
                MPI_Isend(send_buffers_ng[peer_rank].data(), send_buffers_ng[peer_rank].size(), MPI_LONG_LONG,
                        peer_rank, tag, MPI_COMM_WORLD, &req_ng);
                MPI_Isend(send_buffers_ic[peer_rank].data(), send_buffers_ic[peer_rank].size(), MPI_LONG_LONG,
                        peer_rank, tag + 10, MPI_COMM_WORLD, &req_ic);

                requests.push_back(req_ng);
                requests.push_back(req_ic);

                // ========================= Receive =========================
                MPI_Status status_1;
                MPI_Probe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status_1);
                int count;
                MPI_Get_count(&status_1, MPI_LONG_LONG, &count);
                std::vector<int64_t> recvDataID(count);
                MPI_Recv(recvDataID.data(), count, MPI_LONG_LONG, status_1.MPI_SOURCE, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                recv_map_id[static_cast<int64_t>(status_1.MPI_SOURCE)] = recvDataID;

                // ========================= Receive IC Data =========================
                MPI_Status status_2;
                MPI_Probe(MPI_ANY_SOURCE, tag + 10, MPI_COMM_WORLD, &status_2);
                MPI_Get_count(&status_2, MPI_LONG_LONG, &count);
                std::vector<int64_t> recvDataIC(count);
                MPI_Recv(recvDataIC.data(), count, MPI_LONG_LONG, status_2.MPI_SOURCE, tag + 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                recv_map_ic[static_cast<int64_t>(status_2.MPI_SOURCE)] = recvDataIC;
            }

            // ========================= Wait for all sends to complete =========================
            if (!requests.empty())
                MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
        }
    };

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

        if (rank < remainder) 
        {
            localGridSize = cellsPerProc + 1;
            startIdx = rank * localGridSize;
        } else 
        {
            localGridSize = cellsPerProc;
            startIdx = remainder * (cellsPerProc + 1) + (rank - remainder) * cellsPerProc;
        }

        endIdx = startIdx + localGridSize;
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

    std::pair<int64_t, int64_t> globalToLocal(int64_t globalIndex, int64_t gridSize) 
    {
        int numProcs;
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

        if (globalIndex < 0 || globalIndex >= gridSize)
            throw std::out_of_range("globalIndex out of range");

        int64_t cellsPerProc = gridSize / numProcs;
        int64_t remainder    = gridSize % numProcs;

        if (globalIndex < (cellsPerProc + 1) * remainder) 
        {
            int64_t block = cellsPerProc + 1;
            int64_t rank = globalIndex / block;
            int64_t localIndex = globalIndex % block;
            return {rank, localIndex};
        } 
        else 
        {
            int64_t g = globalIndex - remainder * (cellsPerProc + 1);
            int64_t rank = remainder + g / cellsPerProc;
            int64_t localIndex = g % cellsPerProc;
            return {rank, localIndex};
        }
    }

    void getMPICommunicationIDS(const int64_t startID, const int64_t endID, const int64_t gridSize,
                                const std::vector<std::vector<int64_t>> &gridConnect,
                                const std::vector<int> &bndTypes,
                                std::vector<int64_t> &comm_ids, 
                                std::vector<int64_t> &comm_ids_ic,
                                std::vector<int64_t> &comm_ids_ng, 
                                std::vector<int64_t> &comm_rank)
    {
        
        int rank, numProcs;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

        for (int64_t i = 0; i < gridConnect.size(); i++)
            for (int ic = 0; ic < Q; ic++)
                if (gridConnect[i][ic] >= 0 && (gridConnect[i][ic] >= endID || gridConnect[i][ic] < startID))
                {
                    bool isBoundary = (bndTypes[i] != -1); // has at least 1 boundary connection
                    bool isStream = !isBoundary || _GLOBAL_::ifStream[bndTypes[i]][ic];
                    if (isStream)
                    {
                        comm_ids.push_back(i);
                        comm_ids_ic.push_back(ic);
                        auto [peer, ng] = globalToLocal(gridConnect[i][ic], gridSize);
                        comm_ids_ng.push_back(ng);
                        comm_rank.push_back(peer);
                    }
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

    template <typename T>
    std::vector<T> gatherToRoot(const std::vector<T>& localVec)
    {
        int rank, numProcs;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

        int localSize = static_cast<int>(localVec.size());

        std::vector<int> recvCounts, displs;
        if (rank == 0) 
            recvCounts.resize(numProcs);

        MPI_Gather(&localSize, 1, MPI_INT,
                recvCounts.data(), 1, MPI_INT,
                0, MPI_COMM_WORLD);
        
        std::vector<T> globalVec;
        if (rank == 0) 
        {
            displs.resize(numProcs);
            displs[0] = 0;
            for (int i = 1; i < numProcs; i++)
                displs[i] = displs[i - 1] + recvCounts[i - 1];

            globalVec.resize(displs.back() + recvCounts.back());
        }

        MPI_Datatype dtype = MpiType<T>::type();
        MPI_Gatherv(localVec.data(), localSize, dtype,
                    globalVec.data(), recvCounts.data(), displs.data(), dtype,
                    0, MPI_COMM_WORLD);

        return globalVec;
    }
} // namespace mantiis_parallel
