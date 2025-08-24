# /*
#  * Copyright (c) January 2024
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
#  * Proprietary Software - All Rights Reserved
#  *
#  * Unauthorized copying, modification, or distribution of this software,
#  * or any portion of it is prohibited
#  */

CXX := mpic++

CXXFLAGS := -std=c++17 -MMD -MP
LDFLAGS := 

ALGLIB_DIR:=$(HOME)/Programs/alglib/alg/src
OPENCV_DIR=/usr/local
OPENCV_INCLUDE:=-I $(OPENCV_DIR)/include/opencv4 
OPENCV_LIBS_DIR:=-L $(OPENCV_DIR)/lib
OPENCV_LIBS=-lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_imgcodecs -lopencv_ximgproc

OPENMPI_INCLUDE := -I /usr/local/include
OPENMPI_LIBS_DIR := -L /usr/local/lib
OPENMPI_LIBS := -lmpi

BUILD_TYPE := release
CXXFLAGS += -Ofast -fopenmp

ifeq ($(MAKECMDGOALS),debug)
    CXXFLAGS := -std=c++17 -g -O0 -MMD -MP -I$(ALGLIB_DIR)
    BUILD_TYPE := debug
endif

TARGET := mantiis
SRC := mantiis.cpp $(wildcard $(ALGLIB_DIR)/*.cpp)
OBJ := $(SRC:.cpp=.o)
DEP := $(OBJ:.o=.d)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OPENCV_INCLUDE) $(OPENMPI_INCLUDE) $(OBJ) -o $@ $(LDFLAGS) $(OPENCV_LIBS_DIR) $(OPENCV_LIBS) $(OPENMPI_LIBS_DIR) $(OPENMPI_LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(OPENCV_INCLUDE) $(OPENMPI_INCLUDE) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJ) $(DEP)

CXXFLAGS += -I$(ALGLIB_DIR)

debug: $(TARGET)

-include $(DEP)

.PHONY: clean debug

