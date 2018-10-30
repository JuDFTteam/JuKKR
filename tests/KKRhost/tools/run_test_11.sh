#!/usr/bin/env bash

mkdir tests/KKRhost/test_run11_mpi_1_8
cd tests/KKRhost/test_run11_mpi_1_8
ln -s ../test_inputs/test_11_*/* .
rm DTM GMAT
cp -r ../test_inputs/test_11_*/{DTM,GMAT} .
export OMP_NUM_THREADS=1
mpirun -np 8 ../../kkr.x | tee out_kkr
rm -f gmat tmat gref *for* inputcard_generated.txt

cd DTM
pwd && ls && ls -ltr ..
mpirun -np 8 ../../../kkr.x | tee out_kkr
rm -f gmat tmat gref *for* inputcard_generated.txt

cd ../GMAT
mpirun -np 8 ../../../kkr.x | tee out_kkr
rm -f gmat tmat gref *for* inputcard_generated.txt
cd ../


