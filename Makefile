
# LBM
LBM_INCLUDE ?= -I src
OBJS := main.cpp
SRC_DIR := /src

# ALGLIB 
ALGLIB_DIR ?= $(SRC_DIR)/externalLibs/alg
ALGLIB_INCLUDE ?= -I $(ALGLIB_DIR)/src
ALGLIB_LIBS_DIR = -L/$(ALGLIB_DIR)/src
ALGLIB_LIBS ?= -libalglib.a

# OPENCV
OPENCV_DIR ?= /usr/local
OPENCV_INCLUDE ?= -I $(OPENCV_DIR)/incldue/opencv4
OPENCV_LIBS_DIR ?= -L/$(OPENCV_DIR)/lib
OPENCV_LIBS ?= -lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_imgcodecs -lopencv_ximgproc

# C++
CXX := g++
CXXFLAGS := --std=c++17 -Wall -Wextra -Ofast -fopenmp

# Compile and Link
TARGET := main
all: $(TARGET)

# Rule to link the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(OPENCV_LIBS) \
	$(ALGLIB_LIBS) $(ALGLIB_INCLUDE) $(OPENCV_INCLUDE) $(LBM_INCLUDE) -o $(TARGET)