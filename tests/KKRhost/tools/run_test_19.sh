#!/usr/bin/env bash

mkdir tests/KKRhost/test_run19_mpi_2_4
cd tests/KKRhost/test_run19_mpi_2_4
ln -s ../test_inputs/test_19_*/* .
rm bulk
cp -r ../test_inputs/test_19_*/bulk .
cd bulk
export OMP_NUM_THREADS=4
mpirun -np 1 ../../../kkr.x | tee out_kkr
cd ../
rm decifile_bulk_scf; ln -s bulk/decifile decifile_bulk_scf
export OMP_NUM_THREADS=2
mpirun -np 4 ../../kkr.x | tee out_kkr
