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
    Auxiliary Functions
*/

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <unordered_map>

#include "io.hpp"

struct lattice
{
public:
    //                   0  1  2   3  4  5   6   7  8
    std::vector<int> e_x{0, 1, 0, -1, 0, 1, -1, -1, 1};
    std::vector<int> e_y{0, 0, 1, 0, -1, 1, 1, -1, -1};
    // Opposite directions 0  1  2  3  4  5  6  7  8
    std::vector<int> icOpp{0, 3, 4, 1, 2, 7, 8, 5, 6};

public:
    lattice &operator=(const lattice &that)
    {
        if (this != &that)
        {
            e_x = that.e_x;
            e_y = that.e_y;
            icOpp = that.icOpp;
        }
        return *this;
    }
};

struct Index
{
    int row;
    int col;
};

namespace utils
{

    template <typename T>
    int64_t findFirstElement(const std::vector<T> &, const T &);

    // Function to calculate the distance from given point to the nearest wall element in the 2D vector
    /*
    template <typename T>
    double distanceToNearestWall(const std::vector<std::vector<T>> &P, int i, int j, T targetValue, int R = 10)
    {
        double minDistance = std::numeric_limits<double>::infinity();

        for (int row = i - R; row < i + R; ++row)
        {
            for (int col = j - R; col < j + R; ++col)
            {
                if (row >= 0 && row < P.size() && col >= 0 && col <= P[0].size())
                    if (P[row][col] == targetValue)
                    {
                        double distance = sqrt(pow(static_cast<double>(i - row), 2) + pow(static_cast<double>(j - col), 2));

                        minDistance = std::min(minDistance, distance);
                    }
            }
        }

        if (minDistance < 1.0e6)
            return minDistance;
        else
        {
            // std::cout << "Pixel not found icreasing the radius" << std::endl;
            minDistance = distanceToNearestWall(P, i, j, targetValue, 2 * R);
            return minDistance;
        }
    }
    */
    // Function to find the 1D index of the minimum distance from given point to the nearest wall element in the 2D vector
    template <typename T>
    int64_t distanceToNearestWallIndex(const std::vector<std::vector<T>> &P, int i, int j, T targetValue, double D)
    {
        // double minDistance = std::numeric_limits<double>::infinity();
        // std::vector<double> linDist(P.size() * P[0].size(), -1);
        for (int64_t row = i - D - 2; row < i + D + 2; row++)
        {
            for (int64_t col = j - D - 2; col < j + D + 2; col++)
            {
                if (row >= 0 && row < P.size() && col >= 0 && col < P[0].size())
                {
                    double distance = sqrt(pow(static_cast<double>(i - row), 2) + pow(static_cast<double>(j - col), 2));
                    if (abs(distance - D) < 1.0e-4 && P[row][col] == targetValue)
                        return row * P[0].size() + col; // linDist[row * P[0].size() + col] = distance;
                    // minDistance = std::min(minDistance, distance);
                }
            }
        }
        return -1;
    }

    // Function to find the first index of the given value in 1D vector
    template <typename T>
    int64_t findFirstElement(const std::vector<T> &vec, const T &value)
    {
        int64_t index;
        for (int i = 0; i < vec.size(); i++)
        {
            if (vec[i] == value)
            {
                index = static_cast<int>(i);
                break;
            }
        }
        return index;
    }

    // Function to calculate the ratio of the number of "1"s to the number of "0"s in a matrix of integers
    template <typename T> // this function is templated for to int and long compatibility
    double calculateRatio(const std::vector<std::vector<T>> &matrix, T &val)
    {
        int64_t countVal = 0;

        for (int i = 0; i < matrix.size(); i++)
        {
            for (int j = 0; j < matrix[0].size(); j++)
            {
                if (matrix[i][j] == val)
                    countVal++;
            }
        }
        double totalCount = static_cast<double>(matrix.size() * matrix[0].size());
        return static_cast<double>(countVal) / totalCount;
    }

    // Function to calculate the arithmetic mean value of a 2D vector
    template <typename T>
    double calculateMean(const std::vector<std::vector<T>> &matrix)
    {
        T sum = 0.0;
        int64_t count = 0;

        for (const auto &row : matrix)
        {
            for (const auto &element : row)
            {
                sum += static_cast<T>(element);
                count++;
            }
        }
        return static_cast<double>(sum) / static_cast<double>(count);
    }

    // Function to find the minimum value in a 2D vector
    template <typename T>
    T findMinValue(const std::vector<std::vector<T>> &matrix)
    {
        T minValue = std::numeric_limits<T>::max();

        for (const auto &row : matrix)
        {
            for (const auto &element : row)
            {
                if (element < minValue)
                {
                    minValue = element;
                }
            }
        }
        return minValue;
    }

    // Overlaod function to find the minimum value in a 1D vector
    template <typename T>
    T findMinValue(const std::vector<T> &vec)
    {
        T minValue = std::numeric_limits<T>::max();

        for (const auto &element : vec)
        {
            if (element < minValue)
            {
                minValue = element;
            }
        }

        return minValue;
    }

    // Function to find the maximum value in a 1sD vector
    template <typename T>
    T findMaxValue(const std::vector<T> &vec)
    {
        T maxValue = std::numeric_limits<T>::min();

        for (const auto &element : vec)
        {
            if (element > maxValue)
            {
                maxValue = element;
            }
        }

        return maxValue;
    }

    template <typename Key, typename Value>
    bool containsValue(const std::unordered_map<Key, Value>& map, const Value& value)
    {
        for (const auto& pair : map) 
        {   
            if (pair.second == value)
            {
                return true; // Value found in the map
            }
        }
        return false; // Value not found in the map
    }

    template<typename T>
    int closestElement(const std::vector<T>& vec, const T& value)
    {
        if (vec.empty())
            throw std::invalid_argument("Vector is empty");

        // Initialize closest element to the first element of the vector
        T closest = vec[0];
        // Initialize minimum absolute difference to the difference between the first element and the value
        auto minDiff = std::abs(vec[0] - value);
        int closestIndex = 0;
        // Iterate through the vector starting from the second element
        for (int i = 1; i < vec.size(); i++)
        {
            // Calculate the absolute difference between the current element and the value
            auto diff = std::abs(vec[i] - value);
            // If the current element is closer, update the closest element and the minimum difference
            if (diff < minDiff)
            {
                closest = vec[i];
                closestIndex = i;
                minDiff = diff;
            }
        }

        return closestIndex;
    }

    // Function to assign a given value to all elements of a 2D vector
    template <typename T>
    void assignValueToMatrix(std::vector<std::vector<T>> &matrix, const T &value)
    {
        for (auto &row : matrix)
        {
            for (auto &element : row)
            {
                element = value;
            }
        }
    }

    // Function to get unique values of 2D matrix
    template <typename T>
    std::vector<T> flattenAndUnique(const std::vector<std::vector<T>> &matrix)
    {
        // Flatten the 2D vector to a linear vector
        std::vector<T> flattenedVector;
        for (const auto &row : matrix)
        {
            flattenedVector.insert(flattenedVector.end(), row.begin(), row.end());
        }

        // Sort the vector to make unique values adjacent
        std::sort(flattenedVector.begin(), flattenedVector.end());

        // Remove duplicate values
        auto last = std::unique(flattenedVector.begin(), flattenedVector.end());
        flattenedVector.erase(last, flattenedVector.end());

        return flattenedVector;
    }

    // Function to find all 2D indices of an element in a 2D vector
    template <typename T>
    std::vector<std::vector<int>> findAllElements(const std::vector<std::vector<T>> &matrix, const T &value)
    {
        std::vector<std::vector<int>> indices;

        for (int i = 0; i < matrix.size(); i++)
        {
            for (int j = 0; j < matrix[i].size(); j++)
            {
                if (matrix[i][j] == value)
                {
                    indices.push_back({i, j});
                }
            }
        }
        return indices;
    }

    // Overload function to find all 1D indices of an element in a 1D vector
    template <typename T>
    std::vector<int64_t> findAllElements(const std::vector<T> &vec, const T &value)
    {

        // counting
        int64_t N = 0;
        for (int64_t i = 0; i < vec.size(); i++)
        {
            if (vec[i] == value)
            {
                N++;
            }
        }
        //
        std::vector<int64_t> indices(N);
        int64_t count = 0;
        for (int64_t i = 0; i < vec.size(); i++)
        {
            if (vec[i] == value)
            {
                indices[count] = i;
                count++;
            }
        }
        return indices;
    }

    // Function to convert 1D index of an N-dimensional matrix to subscripts
    std::vector<int> ind2sub(const std::vector<int> &size, int64_t index)
    {
        std::vector<int> subscripts(size.size());

        for (int i = 0; i < size.size(); i++)
        {
            int dimSize = size[i];
            subscripts[i] = index % dimSize;
            index /= dimSize;
        }
        return subscripts;
    }

    // Function to convert N-dimensional matrix subscripts to 1D index
    template <typename T>
    int64_t sub2ind(const std::vector<int> &sizes, const std::vector<T> &subscripts)
    {

        // Check if the number of subscripts and sizes match
        if (subscripts.size() != sizes.size())
        {
            std::cerr << "Error: Number of subscripts and sizes do not match." << std::endl;
            return -1; // Error code
        }

        int64_t index = 0;
        int64_t multiplier = 1;
        // switch to column order
        // std::vector<int> subs = {subscripts[1], subscripts[0], subscripts[2]};
        for (int i = 0; i < subscripts.size(); ++i)
        {
            // Check if subscripts are within bounds
            if (subscripts[i] < 0 || subscripts[i] >= sizes[i])
            {   
                std::cerr << "Error: Subscript out of bounds." << std::endl;
                return -1; // Error code
            }

            index += static_cast<int64_t>(subscripts[i]) * multiplier;
            multiplier *= sizes[i];
        }

        return index;
    }

    // Function to display 2D binary image using OpenCV
    template <typename T>
    void displayBinaryImage(std::vector<std::vector<T>> &domain, std::string imName, bool displayIm)
    {
        cv::Mat image(domain.size(), domain[0].size(), CV_8UC1);
        // Convert the binary image vector to an OpenCV Mat
        for (int i = 0; i < domain.size(); i++)
        {
            for (int j = 0; j < domain[0].size(); j++)
            {
                image.at<uchar>(i, j) = domain[i][j] * 255; // Assuming 0 represents black and 1 represents white
            }
        }
        // write image
        cv::imwrite(imName, image);
        if (displayIm)
        {
            cv::namedWindow("Binary Image", cv::WINDOW_NORMAL);
            cv::imshow("Binary Image", image);
            cv::waitKey(0);
        }
    }

    // Function to perform thinning iteration
    // Credits: https://github.com/bsdnoobz/zhang-suen-thinning/blob/master/thinning.cpp
    // This function is not to be called from elsewhere
    namespace
    {
        void thinningIteration(cv::Mat &img, int iter)
        {
            CV_Assert(img.channels() == 1);
            CV_Assert(img.depth() != sizeof(uchar));
            CV_Assert(img.rows > 3 && img.cols > 3);

            cv::Mat marker = cv::Mat::zeros(img.size(), CV_8UC1);

            int nRows = img.rows;
            int nCols = img.cols;

            if (img.isContinuous())
            {
                nCols *= nRows;
                nRows = 1;
            }

            int x, y;
            uchar *pAbove;
            uchar *pCurr;
            uchar *pBelow;
            uchar *nw, *no, *ne; // north (pAbove)
            uchar *we, *me, *ea;
            uchar *sw, *so, *se; // south (pBelow)

            uchar *pDst;

            // initialize row pointers
            pAbove = NULL;
            pCurr = img.ptr<uchar>(0);
            pBelow = img.ptr<uchar>(1);

            for (y = 1; y < img.rows - 1; ++y)
            {
                // shift the rows up by one
                pAbove = pCurr;
                pCurr = pBelow;
                pBelow = img.ptr<uchar>(y + 1);

                pDst = marker.ptr<uchar>(y);

                // initialize col pointers
                no = &(pAbove[0]);
                ne = &(pAbove[1]);
                me = &(pCurr[0]);
                ea = &(pCurr[1]);
                so = &(pBelow[0]);
                se = &(pBelow[1]);

                for (x = 1; x < img.cols - 1; ++x)
                {
                    // shift col pointers left by one (scan left to right)
                    nw = no;
                    no = ne;
                    ne = &(pAbove[x + 1]);
                    we = me;
                    me = ea;
                    ea = &(pCurr[x + 1]);
                    sw = so;
                    so = se;
                    se = &(pBelow[x + 1]);

                    int A = (*no == 0 && *ne == 1) + (*ne == 0 && *ea == 1) +
                            (*ea == 0 && *se == 1) + (*se == 0 && *so == 1) +
                            (*so == 0 && *sw == 1) + (*sw == 0 && *we == 1) +
                            (*we == 0 && *nw == 1) + (*nw == 0 && *no == 1);
                    int B = *no + *ne + *ea + *se + *so + *sw + *we + *nw;
                    int m1 = iter == 0 ? (*no * *ea * *so) : (*no * *ea * *we);
                    int m2 = iter == 0 ? (*ea * *so * *we) : (*no * *so * *we);

                    if (A == 1 && (B >= 2 && B <= 5) && m1 == 0 && m2 == 0)
                        pDst[x] = 1;
                }
            }

            img &= ~marker;
        }
    }

    // Function for thinning the given binary image
    // Credits: https://github.com/bsdnoobz/zhang-suen-thinning/blob/master/thinning.cpp
    void thinning(cv::Mat &src, cv::Mat &dst)
    {
        dst = src.clone();
        dst /= 255; // Convert to binary image

        cv::Mat prev = cv::Mat::zeros(dst.size(), CV_8UC1);
        cv::Mat diff;

        do
        {
            thinningIteration(dst, 0);
            thinningIteration(dst, 1);
            cv::absdiff(dst, prev, diff);
            dst.copyTo(prev);
        } while (cv::countNonZero(diff) > 0);

        dst *= 255;
    }

    // Function to resize a binary image represented as a 2D vector
    template <typename T>
    void resizeBinaryImage(const std::vector<std::vector<T>> &binaryImage, std::vector<std::vector<T>> &resizedBinaryImage, int newWidth, int newHeight)
    {
        // Convert the input binary image to a CV_8U matrix (OpenCV format)
        cv::Mat inputImage(binaryImage.size(), binaryImage[0].size(), CV_8U);
        for (int64_t i = 0; i < binaryImage.size(); i++)
        {
            for (int64_t j = 0; j < binaryImage[0].size(); j++)
            {
                inputImage.at<uchar>(i, j) = binaryImage[i][j] ? 255 : 0;
            }
        }

        // Resize the image using OpenCV
        cv::Mat resizedImage;
        cv::resize(inputImage, resizedImage, cv::Size(newWidth, newHeight));

        // Threshold the resized image to keep it binary
        cv::threshold(resizedImage, resizedImage, 128, 255, cv::THRESH_BINARY);

        // Convert the resized image back to a 2D vector
        resizedBinaryImage.resize(newHeight, std::vector<T>(newWidth));
        for (int i = 0; i < newHeight; ++i)
        {
            for (int j = 0; j < newWidth; ++j)
            {
                resizedBinaryImage[i][j] = resizedImage.at<uchar>(i, j) == 255 ? 1 : 0;
            }
        }
    }

    // Function to solve cubic equation using Cardano's method
    std::vector<double> solveCubic(double a, double b, double c, double d)
    {
        // ax^3 + bx^2 + cx + 1 = 0
        std::vector<double> realRoots;
        // depressed cubic equation
        // t^3 + pt + q = 0
        // where t = x + b / (3a)
        double p = (3 * a * c - b * b) / (3 * a * a);
        double q = (2 * b * b * b - 9 * a * b * c + 27 * a * a * d) / (27 * a * a * a);
        double delta = -4 * p * p * p - 27 * q * q;

        if (delta < 0)
        {
            // Cardano's solution
            double u1 = -q / 2 + sqrt(q * q / 4 + p * p * p / 27);
            double u2 = -q / 2 - sqrt(q * q / 4 + p * p * p / 27);
            double t = pow(u1, 1.0 / 3.0) + pow(u2, 1.0 / 3.0);
            double realRoot = t - b / (3 * a);
            realRoots.push_back(realRoot);
        }
        else if (delta == 0)
        {
            // TBD
        }
        else
        {
            // TBD trigonometric solution
        }
        return realRoots;
    }

    // Specialized 9-point bicubic interpolation
    // Not meant for general pupose use
    /*
    class BicubicInterpolator
    {
    public:
        std::vector<double> X;
        std::vector<double> Y;
        std::vector<double> Z;

        lattice latt;
        // neighbor indices, hardcoded based on lattice configuration
        std::vector<int> xright{1, -1, 5, 0, 8, -1, 2, 4, -1};
        std::vector<int> xleft{3, 0, 6, -1, 7, 2, -1, -1, 4};
        std::vector<int> ypositive{2, 5, -1, 6, 0, -1, -1, 3, 1};
        std::vector<int> ynegative{4, 8, 0, 7, -1, 1, 3, -1, -1};

        // A x alpha = b
        Eigen::MatrixXd A;
        Eigen::VectorXd alpha;
        Eigen::VectorXd b;

    public:
        BicubicInterpolator(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z) : X(x), Y(y), Z(z)
        {
            assert(X.size() == 9);
            assert(X.size() == Y.size());
            assert(X.size() == Z.size());
            A.setZero(16, 16);
            alpha.resize(16);
            alpha.setZero();
            b.resize(16);
            b.setZero();
            int counter = 0;
            int counter2;

            // 9-point
            for (int i = 0; i < X.size(); i++)
            {
                b[i] = Z[i];
                for (int j = 0; j <= 3; j++)
                    for (int k = 0; k <= 3; k++)
                    {
                        A(i, j * 4 + k) = pow(X[i], j) * pow(Y[i], k);
                    }
            }

            // first order derivatives
            counter2 = X.size();
            for (int ic = 0; ic < latt.e_x.size(); ic++)
            {
                if (latt.e_x[ic] == 0)
                {
                    // central difference, hard-coded based on lattice configuration
                    b[counter2] = (Z[xright[ic]] - Z[xleft[ic]]) / 2;
                    for (int j = 1; j <= 3; j++)
                        for (int k = 0; k <= 3; k++)
                        {
                            A(counter2, j * 4 + k) = pow(X[ic], j - 1) * pow(Y[ic], k);
                        }
                    counter2++;
                }
                if (latt.e_y[ic] == 0)
                {
                    b[counter2] = (Z[ypositive[ic]] - Z[ynegative[ic]]) / 2;
                    for (int j = 0; j <= 3; j++)
                        for (int k = 1; k <= 3; k++)
                        {
                            A(counter2, j * 4 + k) = pow(X[ic], j) * pow(Y[ic], k - 1);
                        }
                    counter2++;
                }
            }

            // second-order derivative (mixed xy)
            // hardcoded based on lattice configuration
            b[counter2] = 1.0 / 4.0 * (Z[5] - Z[8] - Z[6] + Z[7]);
            for (int j = 1; j <= 3; j++)
                for (int k = 1; k <= 3; k++)
                {
                    A(counter2, j * 4 + k) = pow(X[0], j - 1) * pow(Y[0], k - 1);
                }
        }

        void solve()
        {
            // Check if A is square and if its number of rows matches b's size
            if (A.rows() != A.cols() || A.rows() != b.size())
            {
                std::cerr << "Error: Matrix A must be square and its size must match vector b's size.\n";
                exit(1);
            }

            // Solve the linear system using Eigen's built-in solver
            this->alpha = A.fullPivLu().solve(b);
        }

        // Interpolate based on single x and y values
        template <typename T>
        double interpolate(T x, T y)
        {
            double result = 0.0;
            for (int j = 0; j <= 3; j++)
            {
                for (int k = 0; k <= 3; k++)
                {
                    result += alpha(j * 4 + k) * pow(static_cast<double>(x), j) * pow(static_cast<double>(y), k);
                }
            }
            return result;
        }

        // Interpolate based on vectors X and Y
        template <typename T>
        std::vector<double> interpolate(const std::vector<T> &X, const std::vector<T> &Y)
        {
            std::vector<double> result(X.size());
            for (size_t i = 0; i < X.size(); i++)
            {
                result[i] = interpolate(X[i], Y[i]);
            }
            return result;
        }
    };
    */

    void load_system_files(std::vector<std::vector<bool>> &ifStream,
                    std::vector<std::vector<int>> &icsr,
                    std::vector<std::vector<double>> &xy_norm,
                    std::vector<int> &bnd_types)
    {
        using namespace IO;
        
        // Get HOME directory
        const char* home = std::getenv("HOME");
        if (!home) {
            throw std::runtime_error("[ERROR]: HOME environment variable not set.");
        }
        std::string basePath = std::string(home) + "/projects/mantiis/include/system_files/";
        
        // ifStream
        int numberOfTypes = ifStream.size();
        int numberOfDirections = 9; // as in D2Q9
        int totalentries = numberOfTypes*numberOfDirections;
        std::vector<double> data(totalentries, 0.0);
        std::string fName = basePath + "ifstream.txt";
        loadDataFromFile(fName, data);
        for (int i = 0; i < numberOfDirections; i++) // row ID, y axis
        {
            for (int j = 0; j < numberOfTypes; j++) // column ID, x axis
            {
                ifStream[j][i] = static_cast<bool>(data[i * numberOfTypes + j]);
            }
        }
        // icsr
        std::fill(data.begin(), data.end(), 0.0);
        fName = basePath + "icsr.txt";
        loadDataFromFile(fName, data);
        for (int i = 0; i < numberOfDirections; i++) // row ID, y axis
        {
            for (int j = 0; j < numberOfTypes; j++) // column ID, x axis
            {
                icsr[j][i] = static_cast<int>(data[i * numberOfTypes + j]) - 1;
            }
        }
        // xy_norm
        numberOfDirections = 2; // as in 2D x and y
        totalentries = numberOfTypes*numberOfDirections;
        data.resize(totalentries);
        data.shrink_to_fit();
        std::fill(data.begin(), data.end(), 0.0);
        fName = basePath + "xy_norm.txt";
        loadDataFromFile(fName, data);
        for (int i = 0; i < numberOfDirections; i++) // row ID, y axis
        {
            for (int j = 0; j < numberOfTypes; j++) // column ID, x axis
            {
                xy_norm[j][i] = static_cast<double>(data[i * numberOfTypes + j]);
            }
        }

        // boundary types in decimals
        data.resize(numberOfTypes);
        data.shrink_to_fit();
        std::fill(data.begin(), data.end(), 0.0);
        fName = basePath + "bnd_types.txt";
        loadDataFromFile(fName, data);
        for (int j = 0; j < numberOfTypes; j++) // column ID, x axis
        {
            bnd_types[j] = static_cast<int>(data[j]);
        }
    };
}
