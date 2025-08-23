# MANTIIS - Multiscale Adaptive Numerical Transport & Inter-molecular Interaction Solver

This project implements a grid refinement method for the Lattice Boltzmann Method (LBM) to simulate high Knudsen transport. The multi-block approach aims to enhance the performance of the LBM simulations.

## External Dependencies

- **OpenCV**: Used for image processing and computer vision tasks. For installation instructions, please refer to the [OpenCV Documentation](https://docs.opencv.org/4.x/d7/d9f/tutorial_linux_install.html). Make sure to build with opencv_contrib.
- **ALGLIB**: A numerical analysis and data processing library. You can download it from the [ALGLIB Website](https://www.alglib.net/download.php).

## Installation

1. **Clone the Repository**

   ```bash
   git clone git@github.com:nnrustamo/multigrid.git

2. **Build ALGLIB**

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

3. **Make sure directories in Makefile are correct**

    Compile in release
    ```bash
    make -j <n_jobs>

4. **Debugging**

    Compile in debug
    ```bash
    make debug

## Execute

    mpirun -np <n> ./mantiis <num_threads> <is_multiblock> <Nx> <input_folder/> <Nt>

 MPI is still under development.

## Cite our work

 Rustamov, N., Mostaghimi, P., Aryana, S. A. On multi-block lattice Boltzmann method for high Knudsen number flows. Advances in Geo-Energy Research, 2025, 16(2): 143-157. https://doi.org/10.46690/ager.2025.05.06
