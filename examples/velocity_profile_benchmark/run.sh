echo "Compiling..."
g++ main.cpp \
`pkg-config --cflags --libs opencv4` \
-I/src/externalLibs/alg/src src/externalLibs/alg/src/*.cpp \
-Ofast -fopenmp -o main
echo "Running..."
./main
