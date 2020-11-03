#!/usr/bin/env bash

cp -rL test_inputs/test_03.2_NOSOC/ test_run03.2_hybrid_1_3
cd test_run03.2_hybrid_1_3

cd NEWSOSOL_NOSOC/
mpirun -np 3 ../../../kkr.x | tee out_kkr
grep -A10000 'Information on renormalization by Lloyds formula' output.000.txt > out_last.txt
cd ..

cd NEWSOSOL_SOCSCL0/
mpirun -np 3 ../../../kkr.x | tee out_kkr
grep -A10000 'Information on renormalization by Lloyds formula' output.000.txt > out_last.txt
cd ../

cd NEWSOSOL_DECOUPLED_SPINS/
mpirun -np 3 ../../../kkr.x | tee out_kkr
grep -A10000 'Information on renormalization by Lloyds formula' output.000.txt > out_last.txt
cd ../../
