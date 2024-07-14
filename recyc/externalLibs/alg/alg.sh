mkdir lib

g++ -shared -fPIC src/alglibinternal.cpp -o lib/libalglibinternal.so
g++ -shared -fPIC src/alglibmisc.cpp -o lib/libalglibmisc.so
g++ -shared -fPIC src/ap.cpp -o lib/libap.so
g++ -shared -fPIC src/dataanalysis.cpp -o lib/libdataanalysis.so
g++ -shared -fPIC src/diffequations.cpp -o lib/libdiffequations.so
g++ -shared -fPIC src/fasttransforms.cpp -o lib/libfasttransforms.so
g++ -shared -fPIC src/integration.cpp -o lib/libintegration.so
g++ -shared -fPIC src/interpolation.cpp -o lib/libinterpolation.so
g++ -shared -fPIC src/kernels_avx2.cpp -o lib/libkernels_avx2.so
g++ -shared -fPIC src/kernels_fma.cpp -o lib/libkernels_fma.so
g++ -shared -fPIC src/kernels_sse2.cpp -o lib/libkernels_sse2.so
g++ -shared -fPIC src/linalg.cpp -o lib/liblinalg.so
g++ -shared -fPIC src/optimization.cpp -o lib/liboptimization.so
g++ -shared -fPIC src/solvers.cpp -o lib/libsolvers.so
g++ -shared -fPIC src/specialfunctions.cpp -o lib/libspecialfunctions.so
g++ -shared -fPIC src/statistics.cpp -o lib/libstatistics.so
