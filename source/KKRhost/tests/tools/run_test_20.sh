#!/usr/bin/env bash
cd tests
cp -r test_inputs/test_20_godfrin/ test_run20_hybrid_1_3
cd test_run20_hybrid_1_3/godfrinON/
export OMP_NUM_THREADS=1
mpirun -np 3 ../../../kkr.x | tee out_kkr
cd ../godfrinOFF/
mpirun -np 3 ../../../kkr.x | tee out_kkr
cd ../../../
