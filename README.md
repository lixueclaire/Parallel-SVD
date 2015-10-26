parallel SVD Using Jacobis Rotations
===================================
  Parallel SVD using jacobis rotations, implemented in OpenMP.

### generate M*N matrix
  M = # of columns <br />
  N = # of Rows <br />
  Matrix must be squared (M=N) <br />
		python randomMatrix.py M N

### serial algorithm
  -t = print out Timing and # of Iterations <br />
  -p = print out Results (U, S, V) <br />
  -d = Generate the Octave files for debug and verify correctness <br />
		g++ SVD.cpp -o svd
		./svd M N -t -d

### parallel algorithm
		mpic++ -fopenmp OMP_SVD.cpp -o omp_svd
		mpiexec -n 2 ./omp_svd M N -t -d

### validate if the result is right
		g++ Validation.cpp -o validation
		./validation

