#!/usr/bin/env bash

# run prep writeout step
export OMP_NUM_THREADS=8
./kkr.x

# run rhoq step, creates our_rhoq.txt which is compared to reference in verify stage
export OMP_NUM_THREADS=1
mpirun -np 8 ./rhoq.x
