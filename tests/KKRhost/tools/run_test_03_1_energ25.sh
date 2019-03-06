#!/usr/bin/env bash

cp -rL test_inputs/test_03.1/ test_run03.1_energ_hybrid_1_25
cd test_run03.1_energ_hybrid_1_25
mv inputcard_mpienerg inputcard
mpirun -np 25 ../../kkr.x | tee out_kkr
cd ../
