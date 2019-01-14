#!/usr/bin/env bash
cp -r test_inputs/test_21_XCs/ test_run21_hybrid_1_3
export OMP_NUM_THREADS=1

cd test_run21_hybrid_1_3/MJW
mpirun -np 3 ../../../kkr.x | tee out_kkr
cd ../vBH
mpirun -np 3 ../../../kkr.x | tee out_kkr
cd ../VWN
mpirun -np 3 ../../../kkr.x | tee out_kkr
cd ../PW91
mpirun -np 3 ../../../kkr.x | tee out_kkr
cd ../PBE
mpirun -np 3 ../../../kkr.x | tee out_kkr
cd ../PBEsol
mpirun -np 3 ../../../kkr.x | tee out_kkr
cd ../../

