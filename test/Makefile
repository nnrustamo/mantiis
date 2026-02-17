# /*
#  * Copyright (c) 2024 Nijat Rustamov
#  *
#  * Author: Nijat Rustamov
#  * Organization: University of Wyoming
#  * Email: nrustamo@uwyo.edu
#  *
#  * Academic Supervisor: Saman Aryana
#  * Email: saryana@uwyo.edu
#  * Organization: University of Wyoming
#  *
#  * This file is a part of Lattice Boltzmann Simulation Software
#  *
#  * Permission is hereby granted, free of charge, to any person obtaining a copy
#  * of this software and associated documentation files (the "Software"), to deal
#  * in the Software without restriction, including without limitation the rights
#  * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  * copies of the Software, and to permit persons to whom the Software is
#  * furnished to do so, subject to the following conditions:
#  *
#  * The above copyright notice and this permission notice shall be included in all
#  * copies or substantial portions of the Software.
#  *
#  * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  * SOFTWARE.
#  */

CXX := /usr/lib64/openmpi/bin/mpic++

CXXFLAGS := -std=c++20 -MMD -MP
LDFLAGS := 

ALGLIB_DIR:=$(HOME)/Programs/alglib/alg/src
OPENCV_DIR=/usr/local
OPENCV_INCLUDE:=-I $(OPENCV_DIR)/include/opencv4 
OPENCV_LIBS_DIR:=-L $(OPENCV_DIR)/lib64
OPENCV_LIBS=-lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_imgcodecs -lopencv_ximgproc

OPENMPI_INCLUDE := -I /usr/include/openmpi-x86_64
OPENMPI_LIBS_DIR := -L /usr/lib64/openmpi/lib
OPENMPI_LIBS := -lmpi

# Include parent directory for mantiis headers
MANTIIS_INCLUDE := -I ..

BUILD_TYPE := release
CXXFLAGS += -Ofast -fopenmp -pthread

ifeq ($(MAKECMDGOALS),debug)
    CXXFLAGS := -std=c++17 -g -O0 -MMD -MP -I$(ALGLIB_DIR)
    BUILD_TYPE := debug
endif

# Targets
TARGETS := transport_test adsorption_test

# ALGLIB sources
ALGLIB_SRC := $(wildcard $(ALGLIB_DIR)/*.cpp)
ALGLIB_OBJ := $(ALGLIB_SRC:.cpp=.o)

# Test sources and objects
TRANSPORT_OBJ := transport_test.o $(ALGLIB_OBJ)
ADSORPTION_OBJ := adsorption_test.o $(ALGLIB_OBJ)

ALL_OBJ := transport_test.o adsorption_test.o $(ALGLIB_OBJ)
DEP := $(ALL_OBJ:.o=.d)

CXXFLAGS += -I$(ALGLIB_DIR)

.PHONY: all clean debug

all: $(TARGETS)

transport_test: $(TRANSPORT_OBJ)
	$(CXX) $(CXXFLAGS) $(MANTIIS_INCLUDE) $(OPENCV_INCLUDE) $(OPENMPI_INCLUDE) $^ -o $@ $(LDFLAGS) $(OPENCV_LIBS_DIR) $(OPENCV_LIBS) $(OPENMPI_LIBS_DIR) $(OPENMPI_LIBS)

adsorption_test: $(ADSORPTION_OBJ)
	$(CXX) $(CXXFLAGS) $(MANTIIS_INCLUDE) $(OPENCV_INCLUDE) $(OPENMPI_INCLUDE) $^ -o $@ $(LDFLAGS) $(OPENCV_LIBS_DIR) $(OPENCV_LIBS) $(OPENMPI_LIBS_DIR) $(OPENMPI_LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(MANTIIS_INCLUDE) $(OPENCV_INCLUDE) $(OPENMPI_INCLUDE) -c $< -o $@

clean:
	rm -f $(TARGETS) $(ALL_OBJ) $(DEP)

debug: $(TARGETS)

-include $(DEP)

