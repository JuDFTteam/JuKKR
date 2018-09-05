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

print 'test coverage input:', test_coverage

# some global settings
modes = ['omp',  'mpi', 'hybrid']  #['serial', 'omp',  'mpi', 'hybrid']
npara_pairs = [[1,1], [1,4], [4,1], [2,2]] # first entry is OMP_NUM_THREADS second number of MPI ranks
global_options = ''
#global_options = 'source /usr/local/bin/compilervars-12.sh intel64; source /usr/local/intel/mkl/bin/mklvars.sh intel64'

test_systems = ['test_run%0.2i'%(i) for i in range(1,20)]

# define masks of test_systes for exgensive (i.e. non-serial) tests
# key is the test_coverage that enters as input via sys.argv command line argument
test_coverages = {1:[0], 2:[1], 3:[2], 4:[3], 5:[4], 6:[5], 7:[6], 8:[7], 12:[11], 14:[13], 15:[14], 16:[15], 17:[16], 18:[17]}

# check for SOC run and change test_coverage automatically
if test_coverage<-1000:
    SOCrun = True
    test_coverage+=1000
else:
    SOCrun = False

# use mpi only if test_coverage option is set to negative value
if test_coverage<0:
    modes = ['hybrid']
    npara_pairs = [[1,3]]
    if test_coverage in [-12]:
        npara_pairs = [[1,8]]
        modes = ['mpi']
    # for FERMIOUT option nranks<=natom is needed
    if test_coverage in [-6]:
        npara_pairs = [[1,3]]
    # for Dirac at the moment nranks==1 is needed
    if test_coverage in [-16]:
        npara_pairs = [[1,1]]
    test_coverage = -test_coverage

print 'settings:'
print 'modes:', modes
print 'para_pairs:', npara_pairs
print 'testcase', array(test_systems)[test_coverages[test_coverage]]
print 'SOCrun:', SOCrun

# loop over all combinations
for mode in modes:
    for npara in npara_pairs:
        for testcase in test_systems:
            # run on parallel execs in parallel
            run_calc = True
            if mode=='serial' and npara[0]*npara[1]>1:
                run_calc = False
            if mode=='omp' and npara[1]>1:
                run_calc = False
            if mode=='mpi' and npara[0]>1:
                run_calc = False

            # check and determine test coverage
            if mode != 'serial':
                if testcase not in array(test_systems)[test_coverages[test_coverage]]:
                    run_calc = False

            # change to SOC extention(".1")
            if SOCrun:
                testcase+='.1'

            # run calculation
            if run_calc:
                path = testcase+'_'+mode+'_'+str(npara[0])+'_'+str(npara[1])
                if path not in os.listdir('.'):
                    job = 'mkdir '+path
                    print job
                    call(job, shell=True)
                    job = 'cd '+path+'; '
                    job+= 'ln -s ../test_inputs/test_%s_*/* .; '%(testcase.replace('test_run',''))
                    if SOCrun:
                        job = job.replace('_*/*', '/*')
                    if global_options != '':
                        job+= global_options+'; '
                    job+= 'export OMP_NUM_THREADS=%i; mpirun -np %i ../../kkr.x | tee out_kkr'%(npara[0], npara[1])
                    print job
                    call(job, shell=True)
                    job = 'cd '+path+'; rm -f gmat tmat gref *for* inputcard_generated.txt'
                    print job
                    call(job, shell=True)
