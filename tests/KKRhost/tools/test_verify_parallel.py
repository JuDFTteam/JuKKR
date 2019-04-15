#!/usr/bin/env python

#use print('message') instead of print 'message' in python 2.7 as well:
from __future__ import print_function

#import pytest
import os, pprint
from kkrparser_functions import parse_kkr_outputfile
from numpy import mean, std, array, loadtxt


class Test_parallel():
    """
    check results of simple scf tests without SOC
    """
    def test_compare_parallel_modes2(self):
        cmplist = ['serial_', 'omp_', 'mpi_', 'hybrid_']
        l_para = [[1,1], [1,4], [2,2], [4,1], [1,7]]
        for ipara, jpara in l_para:
           for irun in cmplist:
               if 'serial' in irun and ipara in [1] and jpara in [1]:
                   cmplist_new = [irun+str(ipara)+'_'+str(jpara)]
               elif 'omp' in irun and jpara in [1] and ipara in l_para:
                   cmplist_new.append(irun+str(ipara)+'_'+str(jpara))
               elif 'mpi' in irun and ipara in [1] and jpara in l_para:
                   cmplist_new.append(irun+str(ipara)+'_'+str(jpara))
               elif 'hybrid' in irun and ipara in l_para and jpara in l_para and ipara*jpara<=4:
                   cmplist_new.append(irun+str(ipara)+'_'+str(jpara))
        cmplist = cmplist_new
        path00 = 'test_run02_'
        cmp_modes(cmplist, path00)

    def test_7_8_mpiatom_mpienerg(self):
        # compare mpiatom and mpienerg parallelisation scheme
        cmplist = ['test_run02_serial_1_1',
                   'test_run07_hybrid_1_3',
		   'test_run08_hybrid_1_3']
        cmp_modes(cmplist, '')

    def test_9_multinode(self):
        # compare mpi and hybrid runs forparallelization across multiple nodes
        cmplist = ['test_run09_mpi_1_32', 
		   'test_run09_hybrid_1_32',
		   'test_run09_hybrid_4_8',
		   'test_run09_hybrid_8_4']
        cmp_modes(cmplist, '')

    def test_SOC_parallel(self):
        # compare mpiatom and mpienerg parallelisation scheme for SOC run
        cmplist = ['test_run02.1_hybrid_1_3',
                   'test_run07.1_hybrid_1_3',
		   'test_run08.1_hybrid_1_3']
        cmp_modes(cmplist, '')

    def test_3_Si_lloyd_rest_parallelization(self):
        # compare parallel modes to check rest parallelization
        cmplist = ['test_run03.1_hybrid_1_3',
                   'test_run03.1_hybrid_1_8',
		   'test_run03.1_hybrid_1_9',
                   'test_run03.1_energ_hybrid_1_25']
        cmp_modes(cmplist, '')

