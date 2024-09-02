# Geometric Multigrid for Lattice Boltzmann Method

This project implements a geometric multigrid method for the Lattice Boltzmann Method (LBM) to simulate high Knudsen transport. The geometric multigrid approach aims to enhance the performance and accuracy of the LBM simulations.

## External Dependencies

- **OpenCV**: Used for image processing and computer vision tasks. For installation instructions, please refer to the [OpenCV Documentation](https://docs.opencv.org/4.x/d7/d9f/tutorial_linux_install.html). Make sure to build with opencv_contrib.
- **ALGLIB**: A numerical analysis and data processing library. You can download it from the [ALGLIB Website](https://www.alglib.net/download.php).

## Installation

1. **Clone the Repository**

   ```bash
   git clone git@github.com:nnrustamo/multigrid.git

3. **Build ALGLIB**

    Copy the following and create a Makefile for ALGLIB
    ```bash
    CXX := g++
    CXXFLAGS := -std=c++17 -Wall -Wextra
    INCLUDE_DIRS := -Ialg/src
    OBJ_DIR := alg/src
    SRCS := alg/src/alglibinternal.cpp \
            alg/src/alglibmisc.cpp \
            alg/src/ap.cpp \
            alg/src/dataanalysis.cpp \
            alg/src/diffequations.cpp \
            alg/src/fasttransforms.cpp \
            alg/src/integration.cpp \
            alg/src/interpolation.cpp \
            alg/src/kernels_avx2.cpp \
            alg/src/kernels_fma.cpp \
            alg/src/kernels_sse2.cpp \
            alg/src/linalg.cpp \
            alg/src/optimization.cpp \
            alg/src/solvers.cpp \
            alg/src/specialfunctions.cpp \
            alg/src/statistics.cpp
    OBJS := $(patsubst alg/src/%.cpp, alg/src/%.o, $(SRCS))
    TARGET_LIB := libalglib.a
    all: $(TARGET_LIB)
    $(OBJ_DIR)/%.o: $(OBJ_DIR)/%.cpp
        $(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) -c $< -o $@
    $(TARGET_LIB): $(OBJS)
        ar rcs $@ $(OBJS)


    ```
    Build ALGLIB
    ```bash
    make -j <n_jobs>

2. **Make sure directories in Makefile are correct**

    Compile
    ```bash
    make -j <n_jobs>


