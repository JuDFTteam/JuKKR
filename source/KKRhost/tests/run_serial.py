#!/usr/bin/env python

from subprocess import call
import sys, os
from numpy import array


if len(sys.argv)>1:
    try:
        test_coverage = int(sys.argv[1])
    except:
        print 'Error! argument can only be an integer but got {}'.format(sys.argv[1])
        sys.exit()
else:
    test_coverage = 0

# some global settings
modes = ['serial']
npara_pairs = [[1,1]] # first entry is OMP_NUM_THREADS second number of MPI ranks
global_options = ''
#global_options = 'source /usr/local/bin/compilervars-12.sh intel64; source /usr/local/intel/mkl/bin/mklvars.sh intel64'

test_systems = ['test_run01', 'test_run02', 'test_run03', 'test_run04', 'test_run05', 'test_run06']

# define masks of test_systes for exgensive (i.e. non-serial) tests
# key is the test_coverage that enters as input via sys.argv command line argument
test_coverages = {1:[0], 2:[1], 3:[2], 4:[3], 5:[4], 6:[5]}

# loop over all combinations
for mode in modes:
    for npara in npara_pairs:
        for testcase in test_systems:
            # run on parallel execs in parallel
            run_calc = True
            if mode=='serial' and npara[0]*npara[1]>1:
                run_calc = False

            # check and determine test coverage
            if testcase not in array(test_systems)[test_coverages[test_coverage]]:
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
                    job+= 'export OMP_NUM_THREADS=%i; mpirun -np %i ../../kkr.x | tee out_kkr'%(npara[0], npara[1])
                    print job
                    call(job, shell=True)
                    job = 'cd '+path+'; rm -f gmat tmat gref *for* inputcard_generated.txt'
                    print job
                    call(job, shell=True)
