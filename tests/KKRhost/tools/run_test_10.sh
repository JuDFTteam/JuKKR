#!/usr/bin/env bash

mkdir tests/KKRhost/test_run10_mpi_1_8
cd tests/KKRhost/test_run10_mpi_1_8
ln -s ../test_inputs/test_10_*/* .
export OMP_NUM_THREADS=1
mpirun -np 8 ../../kkr.x | tee out_kkr
rm -f gmat tmat gref *for* inputcard_generated.txt

# run calculation with impurity wavefunctions
rm inputcard
ln -s imp/* .
export OMP_NUM_THREADS=1
mpirun -np 8 ../../kkr.x | tee out_kkr
rm -f gmat tmat gref *for* inputcard_generated.txt
cd ../

