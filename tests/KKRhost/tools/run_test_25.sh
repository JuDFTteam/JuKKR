#!/usr/bin/env bash
cp -rL test_inputs/test_25_Bconstr/ test_run25_hybrid_1_8
export OMP_NUM_THREADS=1

cd test_run25_hybrid_1_8/
mpirun -np 8 ../../kkr.x | tee out_kkr
cd ../
