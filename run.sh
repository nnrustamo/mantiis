clear
echo "Compiling..."
g++ main.cpp \
`pkg-config --cflags --libs opencv4` \
-I/src/externalLibs/alg/src src/externalLibs/alg/libalglib.a \
-o main -Ofast -fopenmp

echo "Running..."
./main
