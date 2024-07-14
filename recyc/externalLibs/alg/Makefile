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
SRCS := src/alglibinternal.cpp \
        src/alglibmisc.cpp \
        src/ap.cpp \
        src/dataanalysis.cpp \
        src/diffequations.cpp \
        src/fasttransforms.cpp \
        src/integration.cpp \
        src/interpolation.cpp \
        src/kernels_avx2.cpp \
        src/kernels_fma.cpp \
        src/kernels_sse2.cpp \
        src/linalg.cpp \
        src/optimization.cpp \
        src/solvers.cpp \
        src/specialfunctions.cpp \
        src/statistics.cpp

# Object files
OBJS := $(patsubst alg/src/%.cpp, alg/src/%.o, $(SRCS))

# Target static library
TARGET_LIB := libalglib.o

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

