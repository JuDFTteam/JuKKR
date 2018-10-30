#!/usr/bin/env bash
cd tests
cp -rL test_inputs/test_23_dos/ test_run23_hybrid_1_3
export OMP_NUM_THREADS=1

cd test_run23_hybrid_1_3/bulk_SOC
mpirun -np 3 ../../../kkr.x | tee out_kkr
cd ../slab_noSOC
mpirun -np 3 ../../../kkr.x | tee out_kkr
cd ../../../
