#ifort -mkl -O3 -xHost -vec-report test.f90 && ./a.out | tee out

ifort -mkl -O3 -xHost -qopt-report=5 test_dot_prod.f90
export OMP_NUM_THREADS=1; ./a.out | tee out
ifort -mkl=parallel -qopenmp -O3 -xHost test_dot_prod.f90
export OMP_NUM_THREADS=4; ./a.out | tee out2
./plot_times.py

ifort -mkl -O3 -xHost -qopt-report=5 test_matmul.f90
export OMP_NUM_THREADS=1; ./a.out | tee out1
ifort -mkl=parallel -qopenmp -O3 -xHost test_matmul.f90
export OMP_NUM_THREADS=4; ./a.out | tee out12
./plot_times2.py
