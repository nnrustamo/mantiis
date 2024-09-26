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
    Create 2D rectancle and add obstacles
*/

#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "io.hpp"
#include "utils.hpp"

class Shape
{
public:
    int Nx;
    int Ny;
    double resolution; // meters
    std::vector<std::vector<bool>> domain;
    std::vector<std::vector<double>> Kn;
    std::vector<std::vector<double>> LocPore;

public:
    // initiate empty shape
    Shape(int x, int y, double res);

    void addHorizontalBoundary(int h);

    void displayAll();

    void addCircle(int, int, int);

    void addRectangle(int, int, int, int);

    void writeToText(std::string&);

    void calculateProperties(double, double, bool);

    void loadExistingModel(std::string &);

    ~Shape();
};

// initiate empty shape
Shape::Shape(int x, int y, double res) : Nx(x), Ny(y), resolution(res)
{
    // resize vectors and fill with ones
    domain.resize(Ny, std::vector<bool>(Nx, true));
    LocPore.resize(Ny, std::vector<double>(Nx, 0));
    Kn.resize(Ny, std::vector<double>(Nx, 0));

    utils::load_system_files(_GLOBAL_::ifStream, _GLOBAL_::icsr, _GLOBAL_::xy_norm_template, _GLOBAL_::bnd_types);
}

// Add horizontal boundaries
void Shape::addHorizontalBoundary(int h)
{
    // h = 0 -> top boundary @ domain[i=0][j=0:Nx-1]
    // h = Ny-1 -> bottom boundary @ domain[i=Ny-1][j=0:Nx-1]

    for (int j = 0; j < Nx; j++)
        domain[h][j] = 0;
}

void Shape::addCircle(int R, int Centerx, int Centery)
{
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            if (sqrt(pow(static_cast<double>(i - Centery), 2) + pow(static_cast<double>(j - Centerx), 2)) <= static_cast<double>(R))
            {
                domain[i][j] = 0;
            }
        }
    }
}

void Shape::addRectangle(int Xst, int Xe, int Yst, int Ye)
{
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            if (j >= Xst && j <= Xe && i >= Yst && i <= Ye)
            {
                domain[i][j] = 0;
            }
        }
    }
}

// print all, not recommended to use for large domains
void Shape::displayAll()
{
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            std::cout << domain[i][j] << std::endl;
        }
    }
}

// Calculate LocalPoreSize and Knudsen number distribution
// NOTE:: BECAUSE OPENCV SUCKS, THIS FUNCTION IS NOT RECOMMENDED FOR 
// DOMAINS OVER A MILLION CELLS
void Shape::calculateProperties(double resolution, double mfp, bool dispSkeleton = false)
{   
    
    this->resolution = resolution;
    // ================================================================== Prepping
    // Expand shape in both directions to obtain correct skeleton
    int Nx_exp = Nx + 2 * floor(std::max(Nx, Ny)), Ny_exp = Ny + 2 * floor(std::max(Nx, Ny));
    std::vector<std::vector<bool>> domain_exp(Ny_exp, std::vector<bool>(Nx_exp, 0));
    // Place the domain in the middle of expanded shape
    int startX = (Nx_exp - Nx) / 2;
    int startY = (Ny_exp - Ny) / 2;
#pragma omp parallel for default(shared)
    for (int i = -startY; i < Ny + startY; i++)
    {
        for (int j = -startX; j < Nx + startX; j++)
        {
            if (i >= 0 && i < Ny && j >= 0 && j < Nx)
            {
                domain_exp[i + startY][j + startX] = domain[i][j];
            }
        }
    }

// Complete Expansion
#pragma omp parallel for default(shared)
    for (int i = -startY; i < Ny + startY; i++)
    {
        for (int j = -startX; j < Nx + startX; j++)
        {
            if (i >= 1 && i < Ny - 1)
            {
                if (j < 0 || j >= Nx)
                {
                    domain_exp[i + startY][j + startX] = 1;
                }
            }
        }
    }

    // Convert the binary image to a CV_8UC1 Mat
    cv::Mat binaryMat_exp(Ny_exp, Nx_exp, CV_8UC1);
#pragma omp parallel for default(shared)
    for (int i = 0; i < Ny_exp; i++)
    {
        for (int j = 0; j < Nx_exp; j++)
        {
            binaryMat_exp.at<uchar>(i, j) = domain_exp[i][j] ? 255 : 0;
        }
    }

    // Create a copy of the binary image (thinning modifies the input)
    cv::Mat skeleton = binaryMat_exp.clone();
    // Apply thinning (skeletonization) using Zhang-Suen algorithm
    // cv::ximgproc::thinning(skeleton, skeleton);
    utils::thinning(skeleton, skeleton);

    // Convert the skeleton image to a binary vector and cut the expaned regions - BW
    std::vector<std::vector<bool>> skeletonVector(Ny, std::vector<bool>(Nx));
#pragma omp parallel for default(shared)
    for (int i = -startY; i < Ny + startY; i++)
    {
        for (int j = -startX; j < Nx + startX; j++)
        {
            if (i >= 0 && i < Ny && j >= 0 && j < Nx)
            {
                skeletonVector[i][j] = (skeleton.at<uchar>(i + startY, j + startX) == 255);
            }
        }
    }
    if (dispSkeleton)
    {
        cv::imshow("Skeleton", skeleton);
        cv::waitKey(0);
    }

    // Deallocate
    domain_exp.clear();
    domain_exp.shrink_to_fit();
    binaryMat_exp.release();
    skeleton.release();

    // ================================================================== Distance Transforms

    //  Distance to the wall
    cv::Mat binaryMat(Ny, Nx, CV_8UC1);
#pragma omp parallel for default(shared)
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            binaryMat.at<uchar>(i, j) = domain[i][j] ? 255 : 0;
        }
    }
    std::vector<std::vector<double>> distanceVector(Ny, std::vector<double>(Nx, 0));
    cv::Mat BW;
    cv::distanceTransform(binaryMat, BW, cv::DIST_L2, cv::DIST_MASK_PRECISE);
#pragma omp parallel for default(shared)
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {

            distanceVector[i][j] = BW.at<float>(i, j);
        }
    }

    // Distance to the midaxis
    cv::Mat binSkel(Ny, Nx, CV_8UC1);
#pragma omp parallel for default(shared)
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            binSkel.at<uchar>(i, j) = skeletonVector[i][j] ? 0 : 255;
        }
    }
    cv::Mat BWskel;
    cv::distanceTransform(binSkel, BWskel, cv::DIST_L2, cv::DIST_MASK_PRECISE);
    std::vector<std::vector<double>> distanceToMidaxis(Ny, std::vector<double>(Nx, 0));
#pragma omp parallel for default(shared)

    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {

            distanceToMidaxis[i][j] = BWskel.at<float>(i, j);
        }
    }

    //
    std::vector<std::vector<double>> midAxisDistance(Ny, std::vector<double>(Nx, 0));
#pragma omp parallel for default(shared)
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            if (domain[i][j] == 1)
            {
                midAxisDistance[i][j] = static_cast<double>(distanceVector[i][j]) * static_cast<double>(skeletonVector[i][j]);
            }
        }
    }

    std::vector<std::vector<int64_t>> distanceIndex(Ny, std::vector<int64_t>(Nx, 0));
#pragma omp parallel for default(shared)
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            if (domain[i][j] == 1)
            {
                distanceIndex[i][j] = utils::distanceToNearestWallIndex(skeletonVector, i, j, true, distanceToMidaxis[i][j]);
                // TODO
                if (distanceIndex[i][j] == -1)
                {
                    std::cout<<"distanceToMidaxis[i][j] = "<<distanceToMidaxis[i][j]<<std::endl;
                }
                // ENDTODO
            }
        }
    }

    // ================================================================== Local pore size
    int64_t linInd;
#pragma omp parallel for default(shared) private(linInd)
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            if (domain[i][j] == 1)
            {
                linInd = distanceIndex[i][j];
                LocPore[i][j] = 2.0 * midAxisDistance[linInd / Nx][linInd % Nx];
                Kn[i][j] = mfp / LocPore[i][j] / resolution;
                // for special validation case only.
                // Kn[i][j] = 0.15;
                // LocPore[i][j] = _GLOBAL_::mfp / Kn[i][j]/this->resolution;
            }
        }
    }
}

void Shape::loadExistingModel(std::string &f)
{   
    std::cout<<"Loading input from saved data\n";
    std::vector<double> data(Ny * Nx);

    // domain
    std::string fName = f;
    fName.append("pore.txt");
    IO::loadDataFromFile(fName, data);
    
    for (int i = 0; i < Ny; i++) // row ID, y axis
    {
        for (int j = 0; j < Nx; j++) // column ID, x axis
        {
            domain[i][j] = static_cast<bool>(data[j * Nx + i]);
        }
    }

    // kn
    fName = f;
    fName.append("Kn.txt");
    std::fill(data.begin(), data.end(), 0.0);
    IO::loadDataFromFile(fName, data);

    for (int i = 0; i < Ny; i++) // row ID, y axis
    {
        for (int j = 0; j < Nx; j++) // column ID, x axis
        {
            Kn[i][j] = static_cast<double>(data[j * Nx + i]);
        }
    }

    // localporesize
    fName = f;
    fName.append("localporesize.txt");
    std::fill(data.begin(), data.end(), 0.0);
    IO::loadDataFromFile(fName, data);
    for (int i = 0; i < Ny; i++) // row ID, y axis
    {
        for (int j = 0; j < Nx; j++) // column ID, x axis
        {
            LocPore[i][j] = static_cast<double>(data[j * Nx + i]);
        }
    }
}

void Shape::writeToText(std::string& folder)
{

    std::string fName;
    fName = folder + "dimensions.txt";
    // File Dimensions
    std::vector<int> domainInfo = {Ny, Nx};
    IO::writeVectorToFile(fName, domainInfo);

    // Domain
    fName = folder + "domain.txt";
    std::vector<bool> linDomain(Nx * Ny);
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            linDomain[i * Nx + j] = domain[i][j];
        }
    }
    IO::writeVectorToFile(fName, linDomain);

    // Pore Size
    fName = folder + "locpore.txt";
    std::vector<double> linLocPore(Nx * Ny);
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            linLocPore[i * Nx + j] = LocPore[i][j];
        }
    }
    IO::writeVectorToFile(fName, linLocPore);

    // Kn
    fName = folder + "kn.txt";
    std::vector<double> linKn(Nx * Ny);
    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            linKn[i * Nx + j] = Kn[i][j];
        }
    }
    IO::writeVectorToFile(fName, linKn);
}

// Destructor
Shape::~Shape()
{
}
