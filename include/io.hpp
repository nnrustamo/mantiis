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
            std::cout<<"Tried loading "<<fileName<<std::endl;
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
        std::ofstream file(fileName);

        if (!file.is_open())
        {
            std::cout<<"Tried writing "<<fileName<<std::endl;
            throw std::runtime_error("[ERROR]: Unable to open file for writing.");
        }

        for (auto k = data.begin(); k != data.end(); k++)
        {
            file << std::setprecision(15) << *k << '\n';
        }

        file.close();
    }
}

