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


CXX := g++
CXXFLAGS := -std=c++17 -Ofast -fopenmp
LDFLAGS := 

ALGLIB_DIR:=$(HOME)/Programs/alglib/alg/src
OPENCV_DIR=/usr/local
OPENCV_INCLUDE:=-I $(OPENCV_DIR)/include/opencv4 
OPENCV_LIBS_DIR:=-L $(OPENCV_DIR)/lib
OPENCV_LIBS=-lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_imgcodecs -lopencv_ximgproc

# Default to release build
BUILD_TYPE := release

# If debug is specified, add debug flags
ifeq ($(MAKECMDGOALS),debug)
    CXXFLAGS := -std=c++17 -g -O0 -fopenmp
    BUILD_TYPE := debug
endif

TARGET := main
SRC := main.cpp $(wildcard $(ALGLIB_DIR)/*.cpp)
OBJ := $(filter-out main.o, $(SRC:.cpp=.o))

$(TARGET): $(OBJ) main.cpp
	@$(CXX) $(CXXFLAGS) $(OPENCV_INCLUDE) main.cpp $(OBJ) -o $@ $(LDFLAGS) $(OPENCV_LIBS_DIR) $(OPENCV_LIBS)

%.o: %.cpp
	@$(CXX) $(CXXFLAGS) $(OPENCV_INCLUDE) -c $< -o $@

clean:
	@rm -f $(TARGET)

CXXFLAGS += -I$(ALGLIB_DIR)

# Add a debug target
debug: $(TARGET)
	@echo "Compiled in debug mode."

.PHONY: clean debug
