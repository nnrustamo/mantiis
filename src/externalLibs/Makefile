# Makefile for compiling individual source files and creating a static library

# Compiler settings
CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra
INCLUDE_DIRS := -Ialg/src


# OpenCV Directory

# ALGLIB Directory


# Directory for object files
OBJ_DIR := alg/src

# Source files
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

# Object files
OBJS := $(patsubst alg/src/%.cpp, alg/src/%.o, $(SRCS))

# Target static library
TARGET_LIB := libalglib.a

# Default target
all: $(TARGET_LIB)

# Rule to compile source files to object files
$(OBJ_DIR)/%.o: $(OBJ_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) -c $< -o $@

# Rule to create static library
$(TARGET_LIB): $(OBJS)
	ar rcs $@ $(OBJS)

# Clean target
clean:
	rm -f $(TARGET_LIB) $(OBJS)
