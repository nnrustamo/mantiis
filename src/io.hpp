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
    File Descrition goes here
*/

#pragma once


#include <vector>
#include <string>
#include <fstream>
#include <iomanip>

namespace IO
{
    /*template <typename T>
    void loadDataFromFile(const std::string &fileName, std::vector<T> &data);*/

    template <typename T>
    void writeVectorToFile(const std::string &fileName, const std::vector<T> &data);

    template <typename T>
    void loadDataFromFile(const std::string &fileName, std::vector<T> &data)
    {
        std::ifstream file(fileName);

        if (!file.is_open())
        {
            throw std::runtime_error("[ERROR]: File not found.");
        }

        for (int64_t i = 0; i < data.size(); i++)
        {
            if (!(file >> data[i]))
            {
                throw std::runtime_error("[ERROR]: Insufficient data in the file.");
            }
        }

        file.close();
    }

    template <typename T>
    void writeVectorToFile(const std::string &fileName, const std::vector<T> &data)
    {
        // Open the file for writing
        std::ofstream file(fileName);

        if (!file.is_open())
        {
            // File couldn't be opened, handle the error as needed
            throw std::runtime_error("[ERROR]: Unable to open file for writing.");
        }

        // Write data to the file
        for (auto k = data.begin(); k != data.end(); k++)
        {
            file << std::setprecision(15) << *k << '\n';
        }

        // Close the file
        file.close();
    }
}

