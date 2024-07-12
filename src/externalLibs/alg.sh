g++ -c alg/src/alglibinternal.cpp -o alg/src/alglibinternal.o
g++ -c alg/src/alglibmisc.cpp -o alg/src/alglibmisc.o
g++ -c alg/src/ap.cpp -o alg/src/ap.o
g++ -c alg/src/dataanalysis.cpp -o alg/src/dataanalysis.o
g++ -c alg/src/diffequations.cpp -o alg/src/diffequations.o
g++ -c alg/src/fasttransforms.cpp -o alg/src/fasttransforms.o
g++ -c alg/src/integration.cpp -o alg/src/integration.o
g++ -c alg/src/interpolation.cpp -o alg/src/interpolation.o
g++ -c alg/src/kernels_avx2.cpp -o alg/src/kernels_avx2.o
g++ -c alg/src/kernels_fma.cpp -o alg/src/kernels_fma.o
g++ -c alg/src/kernels_sse2.cpp -o alg/src/kernels_sse2.o
g++ -c alg/src/linalg.cpp -o alg/src/linalg.o
g++ -c alg/src/optimization.cpp -o alg/src/optimization.o
g++ -c alg/src/solvers.cpp -o alg/src/solvers.o
g++ -c alg/src/specialfunctions.cpp -o alg/src/specialfunctions.o
g++ -c alg/src/statistics.cpp -o alg/src/statistics.o


ar rcs libalglib.a alg/src/alglibinternal.o alg/src/alglibmisc.o alg/src/ap.o alg/src/dataanalysis.o alg/src/diffequations.o alg/src/fasttransforms.o alg/src/integration.o alg/src/interpolation.o alg/src/kernels_avx2.o alg/src/kernels_fma.o alg/src/kernels_sse2.o alg/src/linalg.o alg/src/optimization.o alg/src/solvers.o alg/src/specialfunctions.o alg/src/statistics.o