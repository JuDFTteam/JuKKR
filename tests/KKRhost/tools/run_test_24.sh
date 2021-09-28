#!/usr/bin/env bash
cp -rL test_inputs/test_24_BXCSCL/ test_run24_hybrid_1_3
export OMP_NUM_THREADS=1

cd test_run24_hybrid_1_3/
mpirun -np 3 ../../kkr.x | tee out_kkr
cd ../
