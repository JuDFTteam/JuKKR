#!/usr/bin/env bash
cd tests
cp -r test_inputs/test_22_LDAU/ test_run22_hybrid_1_3
export OMP_NUM_THREADS=1

cd test_run22_hybrid_1_3/noSOC
mpirun -np 3 ../../../kkr.x | tee out_kkr
cd ../SOC
mpirun -np 3 ../../../kkr.x | tee out_kkr
cd ../../../
