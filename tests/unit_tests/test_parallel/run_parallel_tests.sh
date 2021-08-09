#!/usr/bin/env bash

set -x

ulimit -s unlimited
export OMP_STACKSIZE=1G
export I_MPI_DEBUG=12

#source compiler-select intel
#export FI_PROVIDER=tcp
#export PATH=/usr/local/impi/bin:${PATH}
# mpi test
echo 'MPI test'
echo '--------'
mpiifort test_mpi.f90
#/usr/local/impi/bin/mpifort test_mpi.f90
./a.out | tee out_mpi1.txt 
mpirun -np 2 ./a.out | tee out_mpi2.txt 
echo

# OpenMP test
echo 'OpenMP test'
echo '-----------'
ifort -qopenmp test_omp.f90
#/usr/local/impi/bin/mpifort -qopenmp test_omp.f90
export OMP_NUM_THREADS=1 && ./a.out | tee out_omp1.txt
export OMP_NUM_THREADS=2 && ./a.out | tee out_omp2.txt
echo

# hybrid test
echo 'hybrid OpenMP/MPI test'
echo '----------------------'
mpiifort -qopenmp test_hybrid.f90
#/usr/local/impi/bin/mpifort -qopenmp test_hybrid.f90
mpirun -np 2 ./a.out | tee out_hybrid.txt

