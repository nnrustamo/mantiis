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
    Process and construct 2D domain for single-grid and multi-grid system
*/

#pragma once

#include <vector>
#include <array>
#include <iostream>
#include <cmath>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <bitset>
#include <algorithm>

#include "utils.hpp"
#include "shape.hpp"
#include "comm.hpp"

// Quad-data structure
template <typename T>
struct Quadrant
{
    std::vector<std::vector<T>> vec1;
    std::vector<std::vector<T>> vec2;
    std::vector<std::vector<T>> vec3;
    std::vector<std::vector<T>> vec4;

    // Counter to keep track of which vector to pop
    int popCounter;

    Quadrant(const std::vector<std::vector<T>> &v1,
             const std::vector<std::vector<T>> &v2,
             const std::vector<std::vector<T>> &v3,
             const std::vector<std::vector<T>> &v4)
        : vec1(v1), vec2(v2), vec3(v3), vec4(v4), popCounter(1) {}

    std::vector<std::vector<T>> Pop()
    {
        std::vector<std::vector<T>> poppedVector;

        switch (popCounter)
        {
            case 1:
                poppedVector = std::move(vec1);
                break;
            case 2:
                poppedVector = std::move(vec2);
                break;
            case 3:
                poppedVector = std::move(vec3);
                break;
            case 4:
                poppedVector = std::move(vec4);
                break;
            default:
                break;
        }

        // Increment the counter, resetting to 1 if it exceeds 4
        popCounter = (popCounter % 4) + 1;
        return poppedVector;
    }
};

// A node in the data tree represented by struct that can contain child nodes of the the same type
template <typename T>
struct QuadTreeNode
{
    int dep;                               // depth
    bool isLeaf;                           //
    std::vector<QuadTreeNode<T>> children; // vector of nodes
    std::vector<std::vector<T>> leafNode;  // if isLeaf then leafNode will contain a matrix, else leave empty

    QuadTreeNode(int val = 0, bool leaf = false, std::vector<QuadTreeNode<T>> childNodes = {}) : dep(val), isLeaf(leaf), children(childNodes) {}

    void copyToLeafNode(std::vector<std::vector<T>> &matrix)
    {
        int ny = matrix.size();
        int nx = matrix[0].size();
        leafNode.resize(ny);
        for (int i = 0; i < matrix.size(); i++)
        {
            leafNode[i].resize(nx);
            for (int j = 0; j < matrix[0].size(); j++)
                leafNode[i][j] = matrix[i][j];
        }
    }
};

// Store data for coarse and fine buffer for interpolation
struct bufferData
{
public:
    std::unordered_map<int, std::vector<std::vector<int64_t>>> data;

public:
    void add(int key, const std::vector<std::vector<int64_t>> &value)
    {
        data[key] = value;
    }

    std::vector<std::vector<int64_t>> &operator[](int key)
    {
        return data[key];
    }

    std::vector<std::vector<int64_t>> get(int key)
    {
        auto it = data.find(key);
        if (it != data.end())
            return it->second;
        else
            return std::vector<std::vector<int64_t>>(); // Return an empty 2D vector if key not found
    }

    bool exists(int key)
    {
        return data.find(key) != data.end();
    }

    void remove(int key)
    {
        auto it = data.find(key);
        if (it != data.end())
            data.erase(it);
        else
            std::cout << "Key not found, can't remove.\n";
    }
};

class Grid2D
{
public:
    std::vector<int64_t> gridID;
    std::vector<int> gridType;
    std::vector<int> gridLevel;
    std::vector<int> gridBoundaryType;
    std::vector<int> gridIsBuffer;

    std::vector<std::vector<int>> gridIJ;
    std::vector<std::vector<int64_t>> gridConnect;
    std::vector<std::vector<int64_t>> gridParentChild;

    std::vector<int64_t> solidID;
    std::vector<double> locpore;
    std::vector<double> kn;
    std::vector<int> bndTypes;

    int ny;
    int nx;
    int maxLevel;
    lattice latt;
    bufferData coarseBuffer;
    bufferData fineBuffer;
    bool isMultigrid;

    int64_t globalGridSize;
    int64_t localGridSize;
    int64_t startID;
    int64_t endID;
    int64_t cellsPerProc;
   
public:
    Grid2D(const Shape &, const lattice &, bool);

    void initialize(const Shape &);

    void processMultigrid(const Shape &,
                          std::vector<std::vector<int64_t>> &);

    void processSinglegrid(const Shape &, std::vector<std::vector<int64_t>> &);

    void rebuildGrid(const std::vector<std::vector<int64_t>>&,
                         std::vector<int64_t> &,
                         std::vector<int> &,
                         std::vector<bool> &);

    /*
    void gridReindex(const std::vector<std::vector<int64_t>> &,
                     std::vector<int64_t> &,
                     std::vector<int> &,
                     std::vector<int> &,
                     std::vector<bool> &,
                     std::vector<int> &);

    void buildGrid(const std::vector<bool> &,
                   const std::vector<int64_t> &,
                   const std::vector<int> &,
                   const std::vector<int> &,
                   const std::vector<int> &);
    */

    void getProperties(const Shape &);

    void checkProperties(const Shape &);

    void buildConnections(const std::vector<int64_t> &);

    void determineParentChild();

    void orderBufferIndex();

    void getBoundaryTypes();

    // Function to apply quadtree decomposition on a 2D vector
    QuadTreeNode<int> quadtreeDecompose(const std::vector<std::vector<int>> &,
                                        const std::vector<std::vector<double>> &,
                                        const std::vector<std::vector<double>> &,
                                        double &, int, int &, int &, double &);

    // Function to traverse quadtree
    std::vector<std::vector<int64_t>> backwardsQuadTree(const QuadTreeNode<int> &,
                                                        std::vector<std::vector<int64_t>>, int64_t &);

    // Function to split given 2D matrix into 4 quadrants
    template <typename T>
    Quadrant<T> splitIntoQuadrants(const std::vector<std::vector<T>> &);

    template <typename T>
    std::vector<std::vector<T>> concatenateMatrices(const std::vector<std::vector<T>> &,
                                                    const std::vector<std::vector<T>> &,
                                                    const std::vector<std::vector<T>> &,
                                                    const std::vector<std::vector<T>> &);
    
    void MPI_distributeGridData();
    void MPI_broadcast(int root);
    void debugPrintVectors() const;
};

Grid2D::Grid2D(const Shape &ShapeObject, const lattice &latt, bool isMultigrid = false)
    : ny(ShapeObject.Ny), nx(ShapeObject.Nx), globalGridSize(ny * nx), maxLevel(0), latt(latt), isMultigrid(isMultigrid)
{
}

void Grid2D::initialize(const Shape &ShapeObject)
{
    std::vector<std::vector<int64_t>> reconstructedImage;
    std::cout << "Processing Grid..." << std::endl;
    if (isMultigrid)
        processMultigrid(ShapeObject, reconstructedImage);
    else
        processSinglegrid(ShapeObject, reconstructedImage);

    // TODO
    std::string fName = "reconstructed.dat";
    std::vector<int64_t> linDomain(reconstructedImage.size() * reconstructedImage[0].size());
    for (int i = 0; i < reconstructedImage.size(); i++)
        for (int j = 0; j < reconstructedImage[0].size(); j++)
            linDomain[i * reconstructedImage[0].size() + j] = reconstructedImage[i][j];

    IO::writeVectorToFile(fName, linDomain);

    std::vector<int64_t> linInd;
    std::vector<int> gLvl;
    std::vector<bool> pore;

    std::cout << "Rebuilding Grid..." << std::endl; 
    rebuildGrid(reconstructedImage,linInd, gLvl, pore);

    std::cout << "Calcualting Properties..." << std::endl;
    checkProperties(ShapeObject);
    getProperties(ShapeObject);
    std::cout << "Building Connections..." << std::endl;
    buildConnections(linInd);
    std::cout << "Determining Parent-Child..." << std::endl;
    if (isMultigrid)
        determineParentChild();

    std::cout << "Looking for Buffer..." << std::endl;
    if (isMultigrid)
        orderBufferIndex();

    std::cout<< "Calculating boundary types ... "<<std::endl;
    getBoundaryTypes();
}

void Grid2D::processMultigrid(const Shape &ShapeObject,
                              std::vector<std::vector<int64_t>> &reconstructedImage)
{
    // Parent function to call multigrid functions here
    std::vector<std::vector<int>> pore2D(ShapeObject.Ny, std::vector<int>(ShapeObject.Nx, 0));
    std::vector<std::vector<double>> Kn2D(ShapeObject.Ny, std::vector<double>(ShapeObject.Nx, 0));
    for (int i = 0; i < ShapeObject.Ny; i++)
    {
        for (int j = 0; j < ShapeObject.Nx; j++)
        {
            if (ShapeObject.domain[i][j] == 0)
            {
                solidID.push_back(utils::sub2ind(std::vector<int>{ShapeObject.Ny, ShapeObject.Nx}, std::vector<int>{j, i}));
            }
            pore2D[i][j] = static_cast<int>(ShapeObject.domain[i][j]);
            Kn2D[i][j] = ShapeObject.Kn[i][j];
        }
    }

    // TBD Pad the arrays on each dimension
    /*
    auto padding_size = (std::pow(2, std::ceil(std::log2(NY))) - NY) / 2;

    // Check padding size
    if (padding_size < 1)
    {
        throw std::runtime_error("[ERROR]: padding can not be less than one.");
    }
    else if (typeid(padding_size) != typeid(int))
    {
        throw std::runtime_error("[ERROR]: padding is not an integer.");
    }
    // Apply padding
    addZerosAroundMatrix(pore2D, padding_size);
    addZerosAroundMatrix(Kn2D, padding_size);
    addZerosAroundMatrix(LocPore2D, padding_size);
    */

    // Turn 0 (solids) and 1 (voids) into -1 (solids) and 1 (voids)
    double alpha = 0.99; // solid to void ratio threshold for decompoing matrix
    int SOLIDS = -1;
    int VOIDS = 1;

#pragma omp parallel for default(shared)
    for (int i = 0; i < ShapeObject.Nx; i++)
    {
        for (int j = 0; j < ShapeObject.Ny; j++)
        {
            pore2D[i][j] = (pore2D[i][j] == 0) ? SOLIDS : VOIDS;
        }
    }

    // Calculate distance to the nearest wall
    cv::Mat binaryMat(ShapeObject.Ny, ShapeObject.Nx, CV_8UC1);

#pragma omp parallel for default(shared)
    for (int i = 0; i < ShapeObject.Ny; i++)
    {
        for (int j = 0; j < ShapeObject.Nx; j++)
        {
            binaryMat.at<uchar>(i, j) = ShapeObject.domain[i][j] ? 255 : 0;
        }
    }
    std::vector<std::vector<double>> DWall(ShapeObject.Ny, std::vector<double>(ShapeObject.Nx, 0));
    cv::Mat BW;
    cv::distanceTransform(binaryMat, BW, cv::DIST_L2, cv::DIST_MASK_PRECISE);

#pragma omp parallel for default(shared)
    for (int i = 0; i < ShapeObject.Ny; i++)
    {
        for (int j = 0; j < ShapeObject.Nx; j++)
        {

            DWall[i][j] = BW.at<float>(i, j);
        }
    }

    // Decompose
    QuadTreeNode<int> poreTree;
    double res = ShapeObject.resolution * 1.0e9;
    poreTree = quadtreeDecompose(pore2D, Kn2D, DWall, res, 0, SOLIDS, VOIDS, alpha);

    // Backwards decomposition
    int64_t cellCount = 1;
    reconstructedImage.resize(ShapeObject.Ny, std::vector<int64_t>(ShapeObject.Nx, 0));
    reconstructedImage = backwardsQuadTree(poreTree, reconstructedImage, cellCount);
}

void Grid2D::processSinglegrid(const Shape &ShapeObject,
                               std::vector<std::vector<int64_t>> &reconstructedImage)
{

    reconstructedImage.resize(ShapeObject.Ny, std::vector<int64_t>(ShapeObject.Nx, 0));
    int64_t count = 1;
    for (int i = 0; i < ShapeObject.Ny; i++)
    {
        for (int j = 0; j < ShapeObject.Nx; j++)
        {
            if (ShapeObject.domain[i][j] == 0)
            {
                reconstructedImage[i][j] = (-1) * count;
                solidID.push_back(utils::sub2ind(std::vector<int>{ShapeObject.Ny, ShapeObject.Nx}, std::vector<int>{j, i}));
            }
            else
                reconstructedImage[i][j] = count;

            count++;
        }
    }

    // Solid nodes
}


void Grid2D::rebuildGrid(const std::vector<std::vector<int64_t>> &matrix,
                         std::vector<int64_t> &linInd,
                         std::vector<int> &gLvl,
                         std::vector<bool> &pore)
{   
    // ------------------------------------------------------------------------- Linearize indices
    std::cout << "Looking for Maximum Grid Level... " << std::endl;
    // unique elements
    std::vector<int64_t> uniquElements;
    uniquElements = utils::flattenAndUnique(matrix);
    int64_t numCells = uniquElements.size();

    // Returning data
    linInd.resize(numCells);
    std::fill(linInd.begin(), linInd.end(), 0);
    gLvl.resize(numCells);
    std::fill(gLvl.begin(), gLvl.end(), 0);
    pore.resize(numCells);
    std::fill(pore.begin(), pore.end(), 0);
    
    // Intermediate variables
    int64_t cellValue = 0;
    int celLevel = 0;

    // Convert matrix to unordered maps
    std::unordered_map<int64_t, std::vector<Index>> indices;
    int numRows = matrix.size();
    int numCols = matrix[0].size();

    for (int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < numCols; ++j)
        {
            indices[matrix[i][j]].push_back({i, j});
        }
    }

    // --- Max Level
#pragma omp parallel for default(shared) private(cellValue, celLevel)
    for (int i = 0; i < numCells; i++)
    {
        cellValue = uniquElements[i];
        std::vector<Index> IndexOfcellValue = indices[cellValue];
        celLevel = log2(sqrt(IndexOfcellValue.size())) + 1;
        if (celLevel > maxLevel)
            maxLevel = celLevel;
    }

    // --- Linearize Indices
    /*
        cell/grid dimensions/levels are counted from one
    */
    int64_t count = 0;
    std::vector<int64_t> cellsCounted;
    std::vector<int> vecSize = {ny, nx};
    std::vector<int> newVecSize = {ny, nx, maxLevel}; // this vector will be used later in this scope
    
    std::unordered_map<int64_t, std::vector<int>> ij_to_index_map;
    
    int row = 0, col = 0;
    cellValue = 0; celLevel = 0;
    std::cout << "Linearizing Indices... " << std::endl;
// pragma omp parallel for default(shared) private(cellValue, celLevel)
    for (int i = 0; i < ny; i++)
    {
        std::vector<Index> cellIndices = indices[0];
        std::vector<int> subs = {0, 0, 0};
        for (int j = 0; j < nx; j++)
        {
            cellValue = matrix[i][j];
            if (maxLevel > 1)
            {
                // if the cell value has not been considered before
                if (std::count(cellsCounted.begin(), cellsCounted.end(), cellValue) == 0)
                {
                    cellIndices = indices[cellValue];
                    celLevel = log2(sqrt(cellIndices.size())) + 1;
                    std::vector<int> Rows(cellIndices.size());
                    std::vector<int> Columns(cellIndices.size());
                    // convert to 2d index
                    for (int iX = 0; iX < cellIndices.size(); iX++)
                    {
                        Rows[iX] = cellIndices[iX].row;
                        Columns[iX] = cellIndices[iX].col;
                    }
                    // linear index
                    if (celLevel > 1)
                    {
                        row = utils::findMinValue(Rows);
                        col = utils::findMinValue(Columns);
                    }
                    else
                    {
                        row = Rows[0];
                        col = Columns[0];
                    }
                    subs = {col, row, celLevel - 1};
#pragma omp critical
                    { 
                        // float shift  = pow(2, celLevel - 2) - 0.5; // coordinate centralizer
                        cellsCounted.push_back(cellValue);
                        gLvl[count] = celLevel;
                        pore[count] = (cellValue < 0) ? 0 : 1;
                        linInd[count] = utils::sub2ind(newVecSize, subs);
                        ij_to_index_map[linInd[count]] = {col, row, celLevel};
                        count++;
                    }
                }
            }
            else
            {
                cellIndices = indices[cellValue];
                celLevel = 1;
                row = cellIndices[0].row;
                col = cellIndices[0].col;
                subs = {col, row, celLevel - 1};
#pragma omp critical
                {   
                    // float shift  = pow(2, celLevel - 2) - 0.5; // coordinate centralizer
                    gLvl[count] = celLevel;
                    pore[count] = (cellValue < 0) ? 0 : 1;
                    linInd[count] = utils::sub2ind(newVecSize, subs);
                    ij_to_index_map[linInd[count]] = {col, row, celLevel};        
                    count++;
                }
            }
        }
    }

    /* Sort Linearized idices*/
    /*
        std::unordered_map<int64_t, int64_t> A;
        A.reserve(numCells);
    #pragma omp parallel for default(shared)
        for (int64_t i = 0; i < numCells; i++)
        {
            A[linInd[i]] = i;
        }
        std::sort(linInd.begin(), linInd.end());
        std::vector<int> gLvl_copy(numCells);
        std::vector<bool> pore_copy(numCells);
        
    #pragma omp parallel for default(shared)
        for (int64_t i = 0; i < numCells; i++)
        {
            int64_t originalIndex = A[linInd[i]];
            gLvl_copy[i] = gLvl[originalIndex];
            pore_copy[i] = pore[originalIndex];
            
        }
        gLvl = gLvl_copy;
        pore = pore_copy;

        gLvl_copy.clear();
        gLvl_copy.shrink_to_fit();
        pore_copy.clear();
        pore_copy.shrink_to_fit();
    */

    // clean up
    vecSize.clear();
    vecSize.shrink_to_fit();
    cellsCounted.clear();
    cellsCounted.shrink_to_fit();

    std::cout<<"Mapping the Domain..."<<std::endl;
    // ------------------------------------------------------------------------- Create Distant Neighbor Map
    std::unordered_map<int64_t, std::vector<int64_t>> DistantNeighborMap;
    std::unordered_map<int64_t, std::vector<int64_t>> MissingKeysMap;
    std::unordered_map<int64_t, std::vector<double>> MissingKeysDistanceMap;
    float baseRange = pow(2, maxLevel - 1);
    float range = baseRange; // * sqrt(2);
        
    for(int64_t i = 0; i < numCells; i++)
    {       
        std::vector<int64_t> missing_keys;
        std::vector<double> missing_keys_distance_map;
        float shift = pow(2, gLvl[i] - 2) - 0.5; // coordinate centralizer
        // retrieve subscripts
        if (ij_to_index_map.count(linInd[i]) > 0)
        {   
            int discreteJump = pow(2, gLvl[i] - 1);
            std::vector<int64_t> NeighborMap;
            std::vector<int> sbs = ij_to_index_map[linInd[i]];
            std::vector<float> subscript = {static_cast<float>(sbs[0] + shift), static_cast<float>(sbs[1] + shift)};
            for (float i_c = subscript[0] - baseRange; i_c <= subscript[0] + baseRange; i_c += discreteJump)
            {
                for (float j_c = subscript[1] - baseRange; j_c <= subscript[1] + baseRange; j_c += discreteJump)
                {   
                    std::vector<float> value {i_c, j_c};
                    int64_t idx;
                    if (i_c - shift < 0 || j_c - shift < 0 ||i_c - shift >= nx || j_c - shift >= ny)
                        idx = -1;
                    else
                        idx = utils::sub2ind(newVecSize, std::vector<int>{static_cast<int>(value[0] - shift), static_cast<int>(value[1] - shift), gLvl[i] - 1});
                    if (ij_to_index_map.find(idx) != ij_to_index_map.end())
                    {
                        //  find distance
                        // double distance = sqrt(static_cast<double>(pow(static_cast<double>(value[0] - subscript[0]), 2) + pow(static_cast<double>(value[1] - subscript[1]), 2)));
                        NeighborMap.push_back(static_cast<int64_t>(idx));
                        // NeighborMap.push_back(distance);
                    }
                    else if (value[0] < 0 || value[1] < 0 || value[0] >= nx || value[1] >= ny)
                    {   
                        // (-1) idicates boundary 
                        NeighborMap.push_back(static_cast<int64_t>(-1.0));
                        // NeighborMap.push_back(distance);
                    }
                    else
                    {   
                        double distance = sqrt(static_cast<double>(pow(static_cast<double>(value[0] - subscript[0]), 2) + pow(static_cast<double>(value[1] - subscript[1]), 2)));
                        // missing buffer cell at that level in the grid at that distance
                        NeighborMap.push_back(idx);
                        // NeighborMap.push_back(distance);
                        missing_keys.push_back(idx);
                        missing_keys_distance_map.push_back(distance);
                    }
                }
            }
            DistantNeighborMap[i] = NeighborMap;
            MissingKeysMap[i] = missing_keys;
            MissingKeysDistanceMap[i] = missing_keys_distance_map;
        }
    }

    // ------------------------------------------------------------------------- Determine Buffer
    std::cout<<"Determining Buffer..."<<std::endl;
    //  Construct Distant Neighbor
    std::vector<int64_t> bufferIDS;
    std::vector<std::vector<int>> bufferIJ;
    std::vector<int> bufferLevel;
    std::unordered_map<int64_t, bool> is_buffer;
    std::vector<int64_t> new_ids;
    std::vector<int> new_types;
    int64_t bufferCount = 0;

    for(int64_t i = 0; i < numCells; i++)
    {   
        float shift  = pow(2, gLvl[i] - 2) - 0.5; // coordinate centralizer
        is_buffer[linInd[i]] = false;
        // pull up the map of this neighborhood
        std::vector<int64_t> NeighborMap = DistantNeighborMap[i];
        std::vector<int64_t> missingKeys = MissingKeysMap[i];
        
        std::vector<double> missingKeyDistance = MissingKeysDistanceMap[i];
        bool ifInterface = missingKeys.size() == 0 ? 0: 1;

        // level 2 and higher cells will be scanned 
        if (ifInterface && gLvl[i] > 1)
        {
            // retrieve subscripts
            std::vector<int> sbs = ij_to_index_map[linInd[i]];
            std::vector<float> subscript = {static_cast<float>(sbs[0] + shift), static_cast<float>(sbs[1] + shift)};
            // define range
            float maxI = subscript[0] + range;
            float maxJ = subscript[1] + range;
            float minI = subscript[0] - range;
            float minJ = subscript[1] - range;

            // Retrieve indices of cells within the given range
            std::vector<int64_t> inRangeIndices;
            std::vector<int> inRangeLevels;
            for (int64_t i_c = 0; i_c < numCells; i_c++)
            {   
                const std::vector<int> key_ = ij_to_index_map[linInd[i_c]];
                const std::vector<float> key = {static_cast<float>(key_[0] + shift), static_cast<float>(key_[1] + shift)};
                if (key[0] >= minI && key[0] <= maxI && key[1] >= minJ && key[1] <= maxJ) {
                    inRangeIndices.push_back(linInd[i_c]);
                    inRangeLevels.push_back(key_[2]);
                }
            }
            int min_level = *std::min_element(inRangeLevels.begin(), inRangeLevels.end());
            int max_level = *std::max_element(inRangeLevels.begin(), inRangeLevels.end());

            // Loop thru these cells
            for (int r = 0; r < inRangeIndices.size(); r++)
            {   
                int64_t idx = inRangeIndices[r];
                std::vector<int> sbs_ = ij_to_index_map[idx];
                std::vector<float> ij_ng = {static_cast<float>(sbs_[0] + shift), static_cast<float>(sbs_[1] + shift)};
                
                // break up the higher level cells
                if (std::find(NeighborMap.begin(), NeighborMap.end(), idx) != NeighborMap.end())
                {
                    // shift center
                    float shift_upper  = pow(2, gLvl[i] - 2) - 0.5; 
                    int jump = pow(2, gLvl[i] - 2);
                    // double shift_lower  = pow(2, gLvl[i] - 2 - 1) - 0.5; 
                    // double shift = shift_upper - shift_lower;
                    int i_lower = static_cast<int>(subscript[0] - shift_upper);
                    int j_lower = static_cast<int>(subscript[1] - shift_upper);
                    int i_upper = i_lower + jump; //static_cast<int>(subscript[0] + shift_upper);
                    int j_upper = j_lower + jump; //static_cast<int>(subscript[1] + shift_upper);
                    // create new cells
                    // 1, top left
                    std::vector<int> new_ij = {i_lower, j_lower, gLvl[i] - 2};
                    std::vector<int> new_ij_2d = {i_lower, j_lower};
                    int64_t new_id = utils::sub2ind(newVecSize,  new_ij);
                    if (ij_to_index_map.find(new_id) == ij_to_index_map.end())
                    {
                        ij_to_index_map[new_id] = new_ij_2d;
                        is_buffer[new_id] = true;
                        new_ids.push_back(new_id);
                        new_types.push_back(10*(gLvl[i]-1) + gLvl[i]);
                        bufferCount++;
                    }
                    // 2, top right
                    new_ij = {i_upper, j_lower, gLvl[i] - 2};
                    new_ij_2d = {i_upper, j_lower};
                    new_id = utils::sub2ind(newVecSize,  new_ij);
                    if (ij_to_index_map.find(new_id) == ij_to_index_map.end())
                    {
                        ij_to_index_map[new_id] = new_ij_2d;
                        is_buffer[new_id] = true;
                        new_ids.push_back(new_id);
                        new_types.push_back(10*(gLvl[i]-1) + gLvl[i]);
                        bufferCount++;
                    }
                    
                    // 3, bottom left
                    new_ij = {i_lower, j_upper, gLvl[i] - 2};
                    new_ij_2d = {i_lower, j_upper};
                    new_id = utils::sub2ind(newVecSize,  new_ij);
                    if (ij_to_index_map.find(new_id) == ij_to_index_map.end())
                    {
                       ij_to_index_map[new_id] = new_ij_2d;
                        is_buffer[new_id] = true;
                        new_ids.push_back(new_id);
                        new_types.push_back(10*(gLvl[i]-1) + gLvl[i]);
                        bufferCount++;
                    }
                    // 4, bottom right
                    new_ij = {i_upper, j_upper, gLvl[i] - 2};
                    new_ij_2d = {i_upper, j_upper};
                    new_id = utils::sub2ind(newVecSize,  new_ij);
                    if (ij_to_index_map.find(new_id) == ij_to_index_map.end())
                    {
                        ij_to_index_map[new_id] = new_ij_2d;
                        is_buffer[new_id] = true;
                        new_ids.push_back(new_id);
                        new_types.push_back(10*(gLvl[i]-1) + gLvl[i]);
                        bufferCount++;
                    }
                }
                else if (max_level <= gLvl[i])   // upscale lower level cells
                {
                    int minDistID = 0; 
                    double minDist = std::numeric_limits<double>::max();
                    // search for the closest cell
                    // int closestElement = utils::closestElement(missingKeyDistance, dist);
                    for (int i_m = 0; i_m < missingKeys.size(); i_m++)
                    {
                        std::vector<int> sbs_missing = utils::ind2sub(newVecSize, missingKeys[i_m]);
                        std::vector<float> subs_missing = {static_cast<float>(sbs_missing[0] + shift), static_cast<float>(sbs_missing[1] + shift)};
                        double dist = sqrt(static_cast<double>(pow(ij_ng[0] - subs_missing[0], 2) + pow(ij_ng[1] - subs_missing[1], 2)));
                        if (dist < minDist)
                        {   
                            minDist = dist;
                            minDistID = i_m;
                        }
                    }
                    int64_t new_id = missingKeys[minDistID];
                    // missingKeys.erase(missingKeys.begin() + minDistID); 
                    // add cell
                    std::vector<int> new_ij = utils::ind2sub(newVecSize, new_id);
                    std::vector<int> new_ij_2d = {new_ij[0], new_ij[1]};
                    
                    if (ij_to_index_map.find(new_id) == ij_to_index_map.end() && new_id != -1)
                    {
                        ij_to_index_map[new_id] = new_ij_2d;
                        is_buffer[new_id] = true;
                        new_ids.push_back(new_id);
                        new_types.push_back(10*gLvl[i] + gLvl[i]-1);
                        bufferCount++;
                        missingKeys[minDistID] = -1;
                    }
                }
            }
        }
    }
    std::cout<<"Allocating Memory ...."<<std::endl;
    // ------------------------------------------------------------------------- Allocation and Initialization
    int64_t poreCellCount = std::count(pore.begin(), pore.end(), true);
    this->globalGridSize = poreCellCount + bufferCount;

    this->gridID.reserve(this->globalGridSize);
    this->gridID.resize(this->globalGridSize);
    std::fill(this->gridID.begin(), this->gridID.end(), 0);

    this->gridType.reserve(this->globalGridSize);
    this->gridType.resize(this->globalGridSize);
    std::fill(this->gridType.begin(), this->gridType.end(), 0);
    
    this->gridLevel.reserve(this->globalGridSize);
    this->gridLevel.resize(this->globalGridSize);
    std::fill(this->gridLevel.begin(), this->gridLevel.end(), 0);

    this->gridBoundaryType.reserve(this->globalGridSize);
    this->gridBoundaryType.resize(this->globalGridSize);
    std::fill(this->gridBoundaryType.begin(), this->gridBoundaryType.end(), 0);

    this->gridIsBuffer.reserve(this->globalGridSize);
    this->gridIsBuffer.resize(this->globalGridSize);
    std::fill(this->gridIsBuffer.begin(), this->gridIsBuffer.end(), 0);

    this->gridIJ.reserve(this->globalGridSize);
    this->gridIJ.resize(this->globalGridSize);
    for (int64_t i = 0; i < this->globalGridSize; i++)
    {
        this->gridIJ[i].reserve(2);
        this->gridIJ[i].resize(2);
    }

    this->gridConnect.reserve(this->globalGridSize);
    this->gridConnect.resize(this->globalGridSize);
    for (int64_t i = 0; i < this->globalGridSize; i++)
    {
        this->gridConnect[i].reserve(latt.e_x.size());
        this->gridConnect[i].resize(latt.e_x.size());

    }

    this->locpore.reserve(this->globalGridSize);
    this->locpore.resize(this->globalGridSize);
    std::fill(this->locpore.begin(), this->locpore.end(), 0);

    this->kn.reserve(this->globalGridSize);
    this->kn.resize(this->globalGridSize);
    std::fill(this->kn.begin(), this->kn.end(), 0);

    std::cout<<"Recreating Grid ..."<<std::endl;
    // ------------------------------------------------------------------------- Create grid
    int64_t counter_bufer = 0, counter_cell = 0;
    for (int64_t i = 0; i < numCells; i++)
    {
        // Main cells
        if ( pore[i])
        {
            this->gridID[counter_cell] = linInd[i];
            std::vector<int> subs =  utils::ind2sub(newVecSize, linInd[i]);
            this->gridIJ[counter_cell] = {subs[0], subs[1]};
            this->gridLevel[counter_cell] = gLvl[i];
            this->gridType[counter_cell] = 0;
            this->gridIsBuffer[counter_cell] = 0;
            counter_cell ++;
        }
    }
    for (int64_t i = counter_cell; i < this->globalGridSize; i++)
    {
        // Buffer cells
        if (new_ids.size() > 0)
        {   
            int64_t new_id = new_ids[counter_bufer];
            this->gridID[i] = new_id;
            std::vector<int> subs = utils::ind2sub(newVecSize, new_id);
            this->gridIJ[i] = {subs[0], subs[1]};
            this->gridLevel[i] = subs[2] + 1;
            this->gridType[i] = new_types[counter_bufer];
            this->gridIsBuffer[i] = 1;
            counter_bufer ++;
        }
    }
}

void Grid2D::getProperties(const Shape &ShapeObject)
{

    for (int64_t i = 0; i < globalGridSize; i++)
    {
        locpore[i] = ShapeObject.LocPore[gridIJ[i][1]][gridIJ[i][0]];
        kn[i] = ShapeObject.Kn[gridIJ[i][1]][gridIJ[i][0]];
    }
}

void Grid2D::checkProperties(const Shape & ShapeObject)
{
    for (int64_t i = 0; i< globalGridSize; i++)
    {
        const double H = ShapeObject.LocPore[gridIJ[i][1]][gridIJ[i][0]];
        const double Kn = ShapeObject.Kn[gridIJ[i][1]][gridIJ[i][0]];
        const double critical_dx = 9.56 * std::pow(10.0, -3.0) * 1 / Kn;
        if (ShapeObject.resolution >=critical_dx)
        {
            std::cout<<"[WARNING] Domain dx is larger than critical dx. This may lead to unstable or inacurate  results.\n";
            break;
        }
    }
}

void Grid2D::buildConnections(const std::vector<int64_t> &linInd)
{   
    
    // ============================================================================
    //  Regardless of the type and location, each cell requires 9 conncetions for streaming
    std::vector<int> vecSize = {ny, nx, maxLevel};
    int i, j, i_prev, j_prev, jump, dim;
    int64_t ng;

    // Unordered map for fast index search
    std::unordered_map<int64_t, int64_t> indices;
    for (int64_t c = 0; c < globalGridSize; c++)
    {
        indices[gridID[c]] = c;
    }

    // Convert linInd to unordered_set for faster lookup
    std::unordered_set<int64_t> linIndSet(linInd.begin(), linInd.end());

#pragma omp parallel for default(shared) private(i, j, dim, jump, i_prev, j_prev, ng)
    for (int64_t c = 0; c < globalGridSize; c++)
    {
        i = gridIJ[c][0];
        j = gridIJ[c][1];
        dim = gridLevel[c];
        jump = pow(2, dim - 1);

        for (int ngi = 0; ngi < latt.e_x.size(); ngi++)
        {
            i_prev = i + latt.e_x[ngi] * jump;
            j_prev = j + latt.e_y[ngi] * jump;

            // ensure preiodicity in y direction
            j_prev = (j_prev <= -1) ? j_prev + ny : j_prev;
            j_prev = (j_prev >= ny) ? j_prev - ny : j_prev;

            // ensure preiodicity in x direction
            i_prev = (i_prev <= -1) ? i_prev + nx : i_prev;
            i_prev = (i_prev >= nx) ? i_prev - nx : i_prev;

            // find 1D index
            ng = utils::sub2ind(vecSize, std::vector<int>{i_prev, j_prev, dim - 1});

            // if the ng exists at the same dim
            if (indices.find(ng) != indices.end())
                gridConnect[c][ngi] = indices[ng];
            else // if the ng does not exist
            {
                // it could be a boundary
                if (linIndSet.find(ng) != linIndSet.end())
                    gridConnect[c][ngi] = -1; // marked (-1) indicates boundary hit, boundary treatment will be required
                else
                    gridConnect[c][ngi] = -2; // marked (-2) indicates no neighbor, interpolation will be required
            }
        }
    }
}

void Grid2D::determineParentChild()
{
    // determine the number of cells that require parent-child
    std::vector<int> vecSize = {ny, nx, maxLevel};
    int64_t counter = 0;
    for (int64_t c = 0; c < globalGridSize; c++)
    {
        if (gridLevel[c] > 1)
        {
            counter++;
        }
    }
    gridParentChild.resize(counter);
    for (int64_t i = 0; i < counter; i++)
        gridParentChild[i].resize(5);

    // for every coarse interface cell determine parent-child
    counter = 0;
    int i, j, i_prev, j_prev, jump, dim;
    int64_t c_id, c_ng;
    for (int64_t c = 0; c < globalGridSize; c++)
    {
        if (gridLevel[c] > 1)
        {
            i = gridIJ[c][0];
            j = gridIJ[c][1];
            dim = gridLevel[c];
            jump = pow(2, dim - 1);
            gridParentChild[counter][0] = c;

            // the order of cells are: top left, top right, bottom left, bottom right
            // the order is important for the interpolation
            // top left child
            c_ng = utils::sub2ind(vecSize, std::vector<int>{i, j, dim - 2});
            c_id = -1;
            if (std::count(gridID.begin(), gridID.end(), c_ng) != 0)
                c_id = utils::findFirstElement(gridID, c_ng);
            gridParentChild[counter][1] = c_id;

            // top right child
            c_ng = utils::sub2ind(vecSize, std::vector<int>{i + jump / 2, j, dim - 2});
            c_id = -1;
            if (std::count(gridID.begin(), gridID.end(), c_ng) != 0)
                c_id = utils::findFirstElement(gridID, c_ng);
            gridParentChild[counter][2] = c_id;

            // bottom left child
            c_ng = utils::sub2ind(vecSize, std::vector<int>{i, j + jump / 2, dim - 2});
            c_id = -1;
            if (std::count(gridID.begin(), gridID.end(), c_ng) != 0)
                c_id = utils::findFirstElement(gridID, c_ng);
            gridParentChild[counter][3] = c_id;

            // bottom right child
            c_ng = utils::sub2ind(vecSize, std::vector<int>{i + jump / 2, j + jump / 2, dim - 2});
            c_id = -1;
            if (std::count(gridID.begin(), gridID.end(), c_ng) != 0)
                c_id = utils::findFirstElement(gridID, c_ng);
            gridParentChild[counter][4] = c_id;

            //
            counter++;
        }
    }
}

// Create dictionary style buffer indices
void Grid2D::orderBufferIndex()
{
    // Coarse buffers
    int cell_type;
    for (int i = maxLevel; i >= 2; i--)
    {
        // find cells with the cell type
        cell_type = 10 * i + i - 1;
        int64_t N_coarse = std::count(gridType.begin(), gridType.end(), cell_type);
        std::vector<int64_t> coarse_ids(N_coarse);
        std::vector<std::vector<int64_t>> t_buffer(N_coarse);
        coarse_ids = utils::findAllElements(gridType, cell_type);
        // find child cells
        for (int64_t c = 0; c < coarse_ids.size(); c++)
        {
            int64_t cell_id = coarse_ids[c];
            std::vector<std::vector<int>> indices = utils::findAllElements(gridParentChild, cell_id);
            std::vector<int64_t> child_cell_ids = gridParentChild[indices[0][0]];
            t_buffer[c] = child_cell_ids;
        }
        this->coarseBuffer.add(cell_type, t_buffer);
    }

    // Fine buffers
    for (int i = 1; i < maxLevel; i++)
    {
        // find cells with the cell type
        cell_type = 10 * i + i + 1;
        int64_t N_fine = std::count(gridType.begin(), gridType.end(), cell_type);
        std::vector<int64_t> fine_ids(N_fine);
        std::vector<std::vector<int64_t>> t_buffer(N_fine / 4);
        fine_ids = utils::findAllElements(gridType, cell_type);
        int64_t outerCounter = 0;
        for (int64_t c = 0; c < fine_ids.size(); c += 4)
        {
            std::vector<int64_t> child_cell_ids(4);
            int counter = 0;
            for (int64_t c1 = c; c1 < c + 4; c1++)
            {
                child_cell_ids[counter] = fine_ids[c1];
                counter++;
            }

            // find the parent id
            int64_t parent_id;
            for (int pi = 0; pi < gridParentChild.size(); pi++)
            {
                for (int pj = 1; pj < gridParentChild[pi].size(); pj++)
                {
                    if (gridParentChild[pi][pj] == child_cell_ids[0] ||
                        gridParentChild[pi][pj] == child_cell_ids[1] ||
                        gridParentChild[pi][pj] == child_cell_ids[2] ||
                        gridParentChild[pi][pj] == child_cell_ids[3])
                    {
                        parent_id = gridParentChild[pi][0];
                        goto endofloop;
                    }
                }
            }
        endofloop: // Labeling the end of loop

            child_cell_ids.insert(child_cell_ids.begin(), parent_id);
            t_buffer[outerCounter] = child_cell_ids;
            outerCounter++;
        }
        this->fineBuffer.add(cell_type, t_buffer);
    }
}

// Function to set up bounday types
void Grid2D::getBoundaryTypes()
{
    std::vector<int> d2q9_order_to_new_order{7, 3, 6, 4, 0, 2, 8, 1, 5};
    bndTypes.reserve(globalGridSize);
    bndTypes.resize(globalGridSize);
    std::fill(bndTypes.begin(), bndTypes.end(), -1);

    for (int64_t i = 0; i < globalGridSize; i++)
    {
        if (gridIsBuffer[i] == 0)
        {   
            // assemble binary boundary vector
            std::vector<int64_t> connections = gridConnect[i];
            std::string binaryString;
            for (int ic = 0; ic < d2q9_order_to_new_order.size(); ic++)
            {   
                int conn = connections[d2q9_order_to_new_order[ic]];
                switch (conn)
                {
                case -1:
                    binaryString += '0';
                    break;
                default:
                    binaryString += '1';
                    break;
                }
            }
            // convert to decimalmatlab
            
            std::bitset<32> binary(binaryString);
            int decimalNumber = binary.to_ulong();
            // decimalNumber will be equal to 511 where there is no boundary pixel in the connections
            if (decimalNumber != 511 )
            {
                auto checkFind = std::find(_GLOBAL_::bnd_types.begin(), _GLOBAL_::bnd_types.end(), decimalNumber);
                // check if boundary type is contained
                if (checkFind == _GLOBAL_::bnd_types.end())
                {   
                    throw std::runtime_error("[ERROR]: Boundary type is not contained.");

                }
                else 
                    bndTypes[i] = checkFind - _GLOBAL_::bnd_types.begin();
            }
        }
    }
}

// Function to apply quadtree decomposition on a 2D vector
QuadTreeNode<int> Grid2D::quadtreeDecompose(const std::vector<std::vector<int>> &matrix,
                                            const std::vector<std::vector<double>> &kn,
                                            const std::vector<std::vector<double>> &dwall,
                                            double &dx, int maxDepth,
                                            int &solids, int &voids, double &alpha)

{

    // Split matrix into 4 quadrants
    Quadrant<int> quadrants_matrix = splitIntoQuadrants(matrix);
    Quadrant<double> quadrants_kn = splitIntoQuadrants(kn);
    Quadrant<double> quadrants_dwall = splitIntoQuadrants(dwall);

    // Create a non-leaf node for the current quadrant
    QuadTreeNode<int> rootNode(maxDepth);

    // Variables for deciding if further division is possible
    double beta, cond;
    double dx_prime = 0;
    bool division_flag = false;
    int maxAllowedLevel = 2, lvl;
    int minDist = ceil(pow(2, maxAllowedLevel - 1)*sqrt(2));

    // Loop over each matrix (quadrant)
    for (int i = 0; i < 4; i++)
    {
        division_flag = false;
        std::vector<std::vector<int>> q = quadrants_matrix.Pop();
        std::vector<std::vector<double>> q_kn = quadrants_kn.Pop();
        std::vector<std::vector<double>> q_dwall = quadrants_dwall.Pop();

        // Calculate thresholds
        dx_prime = dx * q.size(); // 9.56 * std::pow(10.0, -3.0) * 1 / calculateMean(q_kn); // correlation
        cond = static_cast<double>(2.0*9.56e-3 * (1.0 / utils::calculateMean(q_kn)));
        beta = utils::calculateRatio(q, solids);

        // Kn based decision
        //division_flag = (dx_prime >= cond) ? true : false;

        // Distance to wall based decision
        if (q.size() > 1 && utils::findMinValue(q_dwall) < minDist * q_dwall.size())
            division_flag = true;
        else if (q.size() == 1)
            division_flag = false;
        /*
        else if (q.size() == 1)
            division_flag = false;
        */
        // Maximum allowed level

        if (log2(q.size()) + 1 > maxAllowedLevel)
            division_flag = true;

        // Divide
        if (((1.0 - alpha) < beta && beta < alpha) || division_flag)
        {
            // continue decomposing the quadrant
            rootNode.children.push_back(quadtreeDecompose(q, q_kn, q_dwall, dx, maxDepth + 1, solids, voids, alpha));
        }
        else if (beta > alpha)
        {
            utils::assignValueToMatrix(q, solids);
            QuadTreeNode<int> childNode(maxDepth + 1);
            childNode.isLeaf = true;
            childNode.copyToLeafNode(q);
            rootNode.children.push_back(childNode);
        }
        else
        {
            utils::assignValueToMatrix(q, voids);
            QuadTreeNode<int> childNode(maxDepth + 1);
            childNode.isLeaf = true;
            childNode.copyToLeafNode(q);
            rootNode.children.push_back(childNode);
        }
    }

    return rootNode;
}

// Function to traverse quad tree
std::vector<std::vector<int64_t>> Grid2D::backwardsQuadTree(const QuadTreeNode<int> &node,
                                                            std::vector<std::vector<int64_t>> matrix, int64_t &cellCount)
{
    if (node.isLeaf)
    {
        // Copy the leaf node data to the corresponding region in the result matrix
        for (int i = 0; i < node.leafNode.size(); i++)
        {
            for (int j = 0; j < node.leafNode[i].size(); j++)
            {
                matrix[i][j] = cellCount * node.leafNode[i][j];
            }
        }
        cellCount++;
        return matrix;
    }
    else
    {
        // Recursively traverse the children and reconstruct the matrix
        Quadrant<int64_t> quadtrants = splitIntoQuadrants(matrix);

        std::vector<std::vector<int64_t>> q1 = backwardsQuadTree(node.children[0], quadtrants.Pop(), cellCount); // top left quadrant
        std::vector<std::vector<int64_t>> q2 = backwardsQuadTree(node.children[1], quadtrants.Pop(), cellCount); // top right quadrant
        std::vector<std::vector<int64_t>> q3 = backwardsQuadTree(node.children[2], quadtrants.Pop(), cellCount); // bottom left quadrant
        std::vector<std::vector<int64_t>> q4 = backwardsQuadTree(node.children[3], quadtrants.Pop(), cellCount); // bottom right quadrant
        std::vector<std::vector<int64_t>> concatQ = concatenateMatrices(q1, q2, q3, q4);
        return concatQ;
    }
}

// Function to split given 2D matrix into 4 quadrants
template <typename T>
Quadrant<T> Grid2D::splitIntoQuadrants(const std::vector<std::vector<T>> &matrix)
{
    // Check if the matrix dimensions are divisible by 2
    if (matrix.size() % 2 != 0 || matrix[0].size() % 2 != 0)
    {
        throw std::invalid_argument("Matrix dimensions must be divisible by 2 in each direction.");
    }

    // Calculate the size of each quadrant
    int quadrantSizeRows = matrix.size() / 2;
    int quadrantSizeCols = matrix[0].size() / 2;

    // Initialize vectors for the quadrants
    std::vector<std::vector<T>> quadrant1(quadrantSizeRows, std::vector<T>(quadrantSizeCols)); // top left
    std::vector<std::vector<T>> quadrant2(quadrantSizeRows, std::vector<T>(quadrantSizeCols)); // top right
    std::vector<std::vector<T>> quadrant3(quadrantSizeRows, std::vector<T>(quadrantSizeCols)); // bottom left
    std::vector<std::vector<T>> quadrant4(quadrantSizeRows, std::vector<T>(quadrantSizeCols)); // bottom right

    // Populate the quadrants
    for (int i = 0; i < quadrantSizeRows; i++)
    {
        for (int j = 0; j < quadrantSizeCols; j++)
        {
            quadrant1[i][j] = matrix[i][j];                                       // top left
            quadrant2[i][j] = matrix[i][j + quadrantSizeCols];                    // top right
            quadrant3[i][j] = matrix[i + quadrantSizeRows][j];                    // bottom left
            quadrant4[i][j] = matrix[i + quadrantSizeRows][j + quadrantSizeCols]; // bottom right
        }
    }

    // Return the quadrants
    return Quadrant<T>(quadrant1, quadrant2, quadrant3, quadrant4);
}

// Function to concatenate 4 matrices together
template <typename T>
std::vector<std::vector<T>> Grid2D::concatenateMatrices(const std::vector<std::vector<T>> &vec1,
                                                        const std::vector<std::vector<T>> &vec2,
                                                        const std::vector<std::vector<T>> &vec3,
                                                        const std::vector<std::vector<T>> &vec4)
{

    // Concatenate matrices
    std::vector<std::vector<T>> result;
    result.resize(vec1.size() * 2, std::vector<int64_t>(vec1[0].size() * 2, 0));

    // top left
    for (int i = 0; i < vec1.size(); i++)
        for (int j = 0; j < vec1[0].size(); j++)
            result[i][j] = vec1[i][j];
    // top right
    for (int i = 0; i < vec2.size(); i++)
        for (int j = 0; j < vec2[0].size(); j++)
            result[i][j + vec2[0].size()] = vec2[i][j];
    // bottom left
    for (int i = 0; i < vec3.size(); i++)
        for (int j = 0; j < vec3[0].size(); j++)
            result[i + vec3.size()][j] = vec3[i][j];
    // bottom right
    for (int i = 0; i < vec4.size(); i++)
        for (int j = 0; j < vec4[0].size(); j++)
            result[i + vec4.size()][j + +vec4[0].size()] = vec4[i][j];

    return result;
}

void Grid2D::MPI_distributeGridData()
{   
    int numProcs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    auto [startIdx, endIdx, local_GridSize, cellsPerProc] = mantiis_parallel::calculateMPIGridPartition(globalGridSize);
    
    this->localGridSize = local_GridSize;
    this->startID= startIdx;
    this->endID = endIdx;
    this->cellsPerProc = cellsPerProc;

    // Distribute 1D vectors
    gridID = mantiis_parallel::getSubDomainVector(gridID, startIdx, endIdx);
    gridType = mantiis_parallel::getSubDomainVector(gridType, startIdx, endIdx);
    gridLevel = mantiis_parallel::getSubDomainVector(gridLevel, startIdx, endIdx);
    gridBoundaryType = mantiis_parallel::getSubDomainVector(gridBoundaryType, startIdx, endIdx);
    gridIsBuffer = mantiis_parallel::getSubDomainVector(gridIsBuffer, startIdx, endIdx);
    locpore = mantiis_parallel::getSubDomainVector(locpore, startIdx, endIdx);
    kn = mantiis_parallel::getSubDomainVector(kn, startIdx, endIdx);
    bndTypes = mantiis_parallel::getSubDomainVector(bndTypes, startIdx, endIdx);

    // Distribute 2D vectors
    gridIJ = mantiis_parallel::getSubDomainVector(gridIJ, startIdx, endIdx);
    gridConnect = mantiis_parallel::getSubDomainVector(gridConnect, startIdx, endIdx);

    std::cout << "[INFO] Rank " << rank << " has " << localGridSize << " cells"<<" out of " << globalGridSize
              << ". Start ID: " << startID << ", End ID: " << endID << std::endl;

}

void Grid2D::MPI_broadcast(int root = 0)
{   
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Bcast(&ny, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&nx, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&maxLevel, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&isMultigrid, 1, MPI_C_BOOL, root, MPI_COMM_WORLD);
    MPI_Bcast(&globalGridSize, 1, MPI_LONG_LONG, root, MPI_COMM_WORLD);

    mantiis_parallel::broadcast_vector(gridID, root);
    mantiis_parallel::broadcast_vector(gridType, root);
    mantiis_parallel::broadcast_vector(gridLevel, root);
    mantiis_parallel::broadcast_vector(gridBoundaryType, root);
    mantiis_parallel::broadcast_vector(gridIsBuffer, root);
    mantiis_parallel::broadcast_vector(gridIJ, root);
    mantiis_parallel::broadcast_vector(gridConnect, root);
    mantiis_parallel::broadcast_vector(gridParentChild, root);
    mantiis_parallel::broadcast_vector(solidID, root);
    mantiis_parallel::broadcast_vector(locpore, root);
    mantiis_parallel::broadcast_vector(kn, root);
    mantiis_parallel::broadcast_vector(bndTypes, root);

    if (isMultigrid && rank != root)
        orderBufferIndex();
}

void Grid2D::debugPrintVectors() const
{   

    int numProcs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    auto print1DVector = [](const auto &vec, const std::string &name) 
    {
        std::cout << name << " (size: " << vec.size() << "): ";
        for (const auto &val : vec)
        {
            std::cout << val << "\n";
        }
        std::cout << std::endl;
    };

    auto print2DVector = [](const auto &vec, const std::string &name) 
    {
        std::cout << name << " (size: " << vec.size() << " x ";
        if (!vec.empty())
            std::cout << vec[0].size();
        else
            std::cout << "0";
        std::cout << "):" << std::endl;
        for (const auto &row : vec)
        {
            for (const auto &val : row)
            {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    };

    for (int proc = 0; proc < numProcs; ++proc)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == proc)
        {
            std::cout << "Debugging Grid2D Vectors (Rank " << rank << "):" << std::endl;
            print1DVector(gridID, "gridID");
            print1DVector(gridType, "gridType");
            print1DVector(gridLevel, "gridLevel");
            print1DVector(gridBoundaryType, "gridBoundaryType");
            print1DVector(gridIsBuffer, "gridIsBuffer");
            print1DVector(locpore, "locpore");
            print1DVector(kn, "kn");
            print1DVector(bndTypes, "bndTypes");
            print1DVector(solidID, "solidID");
            print2DVector(gridIJ, "gridIJ");
            print2DVector(gridConnect, "gridConnect");
            std::cout << std::flush;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}
