#!/bin/sh
export OMP_STACKSIZE=20M
mpirun -np $1 kkr2.exe
