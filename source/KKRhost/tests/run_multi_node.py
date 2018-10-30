#!/usr/bin/env python

from subprocess import call
import sys, os
from numpy import array

# some global settings
modes = ['omp',  'mpi', 'hybrid']  #['serial', 'omp',  'mpi', 'hybrid']
npara_pairs = [[1,32], [4,8], [8,4]] # first entry is OMP_NUM_THREADS second number of MPI ranks
global_options = ''
#global_options = 'source /usr/local/bin/compilervars-12.sh intel64; source /usr/local/intel/mkl/bin/mklvars.sh intel64'

test_systems = ['test_run'+str(i) for i in range(1,9)]



# loop over all combinations
for mode in modes:
    for npara in npara_pairs:
        testcase = 'test_run09'
        # run on parallel execs in parallel
        run_calc = True
        if mode=='serial' and npara[0]*npara[1]>1:
            run_calc = False
        if mode=='omp' and npara[1]>1:
            run_calc = False
        if mode=='mpi' and npara[0]>1:
            run_calc = False

        # run calculation
        if run_calc:
            path = testcase+'_'+mode+'_'+str(npara[0])+'_'+str(npara[1])
            if path not in os.listdir('.'):
                job = 'mkdir '+path
                print job
                call(job, shell=True)
                job = 'cd '+path+'; '
                job+= 'ln -s ../test_inputs/test_%s_*/* .; '%(testcase.replace('test_run',''))
                if global_options != '':
                    job+= global_options+'; '
                job+= 'export OMP_NUM_THREADS=%i; srun --nodes=%i --ntasks-per-node=%i ../../kkr.x | tee out_kkr'%(npara[0], npara[1]/(8/npara[0]), 8/npara[0])
                print job
                call(job, shell=True)
                job = 'cd '+path+'; rm -f gmat tmat gref *for* inputcard_generated.txt'
                print job
                call(job, shell=True)
