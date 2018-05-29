#!/usr/bin/env python

import pytest
import os, pprint
from kkrparser_functions import parse_kkr_outputfile
from numpy import mean, std, array, loadtxt


class Test_check_test_runs():
    """
    check results of the tests
    """

    def test_verify1_Au_bulk(self):
        path0 = 'test_run1_serial_1_1/'
        standard_verify(path0, rms_threshold=5*10**-8)

    def test_verify2_Fe_slab(self):
        path0 = 'test_run2_serial_1_1/'
        standard_verify(path0, rms_threshold=7*10**-8)

    def test_verify3_Si_lloyd(self):
        path0 = 'test_run3_serial_1_1/'
        standard_verify(path0, rms_threshold=8*10**-9)

    """
    def test_compare_parallel_modes1(self):
        cmplist = ['serial_1_1', 'omp_1_1', 'omp_4_1', 'mpi_1_1' ,'mpi_1_4', 'hybrid_1_1', 'hybrid_1_4', 'hybrid_4_1', 'hybrid_2_2']
        path00 = 'test_run1_'
        cmp_modes(cmplist, path00)
    """

    def test_compare_parallel_modes2(self):
        cmplist = ['serial_', 'omp_', 'mpi_', 'hybrid_']
        l_para = [1,2,4,8]
        for ipara in l_para:
            for jpara in l_para:
                for irun in cmplist:
                    if 'serial' in irun and ipara in [1] and jpara in [1]:
                        cmplist_new = [irun+str(ipara)+'_'+str(jpara)]
                    elif 'omp' in irun and jpara in [1] and ipara in l_para:
                        cmplist_new.append(irun+str(ipara)+'_'+str(jpara))
                    elif 'mpi' in irun and ipara in [1] and jpara in l_para:
                        cmplist_new.append(irun+str(ipara)+'_'+str(jpara))
                    elif 'hybrid' in irun and ipara in l_para and jpara in l_para and ipara*jpara<=8:
                        cmplist_new.append(irun+str(ipara)+'_'+str(jpara))
        cmplist = cmplist_new
        path00 = 'test_run2_'
        cmp_modes(cmplist, path00)

    """
    def test_compare_parallel_modes3(self):
        cmplist = ['serial_1_1', 'omp_1_1', 'omp_4_1', 'mpi_1_1' ,'mpi_1_4', 'hybrid_1_1', 'hybrid_1_4', 'hybrid_4_1', 'hybrid_2_2']
        path00 = 'test_run3_'
        cmp_modes(cmplist, path00)
    """

    def test_verify4_Jijs_noSOC(self):
        paths = 'test_run4_serial_1_1/ test_run4_mpi_1_4/ test_run4_hybrid_1_4/'.split()
        path0 = 'test_inputs/test_4_Jijs_Fe_slab_lmax2_noSOC/ref/'
        # compare Jij.atom* files of different runs with reference (in path0)
        for path in paths:
           files = 'Jij.atom00002 Jij.atom00003'.split()
           for fname in files:
              num, text = read_file(path+fname)
              num_ref, text_ref = read_file(path0+fname)
              assert std(num-num_ref)<10**-10
              assert set(text)-set(text_ref)==set()

    def test_verify5_kkrflex(self):
        paths = 'test_run5_serial_1_1/ test_run5_mpi_1_4/ test_run5_hybrid_1_4/'.split()
        path0 = 'test_inputs/test_5_Silicon_lloyd_kkrflex_output_lmax2_noSOC/ref/'
        # compare kkrflex_* files of different runs with reference (in path0)
        for path in paths:
           files = 'kkrflex_atominfo kkrflex_hoststructure.dat kkrflex_intercell_cmoms kkrflex_intercell_ref kkrflex_tmat'.split()
           for fname in files:
              num, text = read_file(path+fname)
              num_ref, text_ref = read_file(path0+fname)
              assert std(num-num_ref)<10**-10
              assert set(text)-set(text_ref)==set()

    def test_verify6_FERMIOUT(self):
        paths = 'test_run6_serial_1_1/ test_run6_mpi_1_4/ test_run6_hybrid_1_4/'.split()
        path0 = 'test_inputs/test_6_Silicon_lloyd_FERMIOUT_output_lmax2_noSOC/ref/'
        # compare TBkkr_* files of different runs with reference (in path0)
        for path in paths:
           files = 'TBkkr_container.txt  TBkkr_params.txt'.split()
           for fname in files:
              num, text = read_file(path+fname)
              num_ref, text_ref = read_file(path0+fname)
              assert std(num-num_ref)<10**-10
              assert set(text)-set(text_ref)==set()

    def test_verify7_8_mpiatom_mpienerg(self):
        # compare mpiatom and mpienerg parallelisation scheme
        cmplist = ['test_run7_mpi_1_2', 'test_run7_hybrid_1_2', 
		   'test_run7_mpi_1_3', 'test_run7_hybrid_1_3',
		   'test_run7_mpi_1_4', 'test_run7_hybrid_1_4',
                   'test_run8_mpi_1_2', 'test_run8_hybrid_1_2', 
		   'test_run8_mpi_1_3', 'test_run8_hybrid_1_3',
		   'test_run8_mpi_1_4', 'test_run8_hybrid_1_4']
        cmplist+= ['test_run7_mpi_1_7', 'test_run7_hybrid_1_7', 
		   'test_run7_mpi_1_8', 'test_run7_hybrid_1_8',
		   'test_run8_mpi_1_7', 'test_run8_hybrid_1_7',
		   'test_run8_mpi_1_8', 'test_run8_hybrid_1_8']
        cmp_modes(cmplist, '')

    def test_verify11_multinode(self):
        # compare mpi and hybrid runs forparallelization across multiple nodes
        cmplist = ['test_run11_mpi_1_32', 
		   'test_run11_hybrid_1_32',
		   'test_run11_hybrid_4_8',
		   'test_run11_hybrid_8_4']
        cmp_modes(cmplist, '')

    def test_verify12_OPERATOR(self):
        path  = 'test_run12_mpi_1_8/'
        path0 = 'test_run12_mpi_1_8/ref/'
        # compare TBkkr_rhod.txt file with reference (in path0)
        fname = 'TBkkr_rhod.txt'
        num, text = read_file(path+fname)
        num_ref, text_ref = read_file(path0+fname)
        assert std(abs(num-num_ref))<10**-14
        assert mean(abs(num-num_ref))<10**-14
        assert abs(num-num_ref).max()<5*10**-13
        assert set(text)-set(text_ref)==set()
        # compare output of OPERATOR for host and for impurity wavefunctions
        for filename in 'TBkkr_rhod.txt TBkkr_torq.txt TBkkr_spinflux.txt'.split():
          d = loadtxt(path+filename)
          d0 = loadtxt(path+filename.replace('.txt', '_imp.txt'))
          nsigma = 3
          if 'rhod' in filename:
             nsigma +=1
          d1 = d[:,0].reshape(nsigma,72, 32, 32); d1 = d1[:,36:48,:,:]; d01 = d0[:,0].reshape(nsigma,12,32,32)
          d2 = d[:,1].reshape(nsigma,72, 32, 32); d2 = d2[:,36:48,:,:]; d02 = d0[:,1].reshape(nsigma,12,32,32)
          d1 = d1.reshape(-1); d2 = d2.reshape(-1); d01 = d01.reshape(-1); d02 = d02.reshape(-1)
          diff1 = d01-d1; diff2 = d02-d2
          assert mean(diff1) < 10**-15
          assert abs(diff1).max() < 10**-15
          assert mean(diff2) < 10**-15
          assert abs(diff2).max() < 10**-15

    def test_verify13_DTM_GMAT(self):
        path  = 'test_run13_mpi_1_8/'
        path0 = 'test_run13_mpi_1_8/ref/'
        for f in 'DTM/DTMTRX ./green_host GMAT/GMATLL_GES'.split():
           fname = f
           num, text = read_file(path+fname)
           num_ref, text_ref = read_file(path0+fname.split('/')[1])
           print fname
           print std(abs(num-num_ref))
           print mean(abs(num-num_ref))
           print abs(num-num_ref).max()
           print set(text)-set(text_ref)==set()
           assert std(abs(num-num_ref))<5*10**-11
           assert mean(abs(num-num_ref))<10**-12
           assert abs(num-num_ref).max()<2*10**-8
           assert set(text)-set(text_ref)==set()

    def test_verify14_qdos(self):
        path  = 'test_run14_mpi_1_8/'
        path0 = 'test_run14_mpi_1_8/ref/'
        for f in 'qdos.01.1.dat qdos.01.2.dat qdos.02.1.dat qdos.02.2.dat qdos.03.1.dat qdos.03.2.dat qdos.04.1.dat qdos.04.2.dat'.split():
           fname = f
           num, text = read_file(path+fname)
           num_ref, text_ref = read_file(path0+fname)
           # remove line with serial number
           text = text[1:]
           text_ref = text_ref[1:]
           # now compare
           print fname
           print std(abs(num-num_ref))
           print mean(abs(num-num_ref))
           print abs(num-num_ref).max()
           print set(text)-set(text_ref)==set()
           assert std(abs(num-num_ref))<5*10**-16
           assert mean(abs(num-num_ref))<10**-14
           assert abs(num-num_ref).max()<2*10**-12
           assert set(text)-set(text_ref)==set()

        
# helper functions

def standard_verify(path0, rms_threshold=10**-8, rms_threshold_end=10**-8, neutr_threshold=10**-6):
    """
    wrapper for standard tests reading output and comparins rms and charge neutrality
    """
    # use parser function from aiida-kkr
    success, parser_msgs, out_dict = parse_kkr_outputfile({}, path0+'out_kkr', path0+'output.0.txt', path0+'output.000.txt', path0+'out_timing.000.txt', path0+'out_potential', path0+'nonco_angle_out.dat')
    pprint.pprint(parser_msgs)
    pprint.pprint(out_dict)
    # first check if parsing was successful
    assert success
    # check if initial iteration is still converged
    assert out_dict['convergence_group']['rms_all_iterations'][0] < rms_threshold
    # check if iteration procedure was successful
    assert out_dict['convergence_group']['rms'] < rms_threshold_end
    # check if charge neutrality is correct
    assert abs(out_dict['convergence_group']['charge_neutrality']) <= neutr_threshold

def cmp_modes(cmplist, path00):
    """
    check convergence and charges across parallel runs
    returns dict with rms and charges entries for all paths given in 'cmplist' in parent directory 'path00'
    """
    cmp_values = {'rms':[], 'charges':[]}
    for addpath in cmplist:
        path0 = path00+addpath
        print path0, os.listdir('.')
        if path0 in os.listdir('.'):
            path0+='/'
            success, parser_msgs, out_dict = parse_kkr_outputfile({}, path0+'out_kkr', path0+'output.0.txt', path0+'output.000.txt', path0+'out_timing.000.txt', path0+'out_potential', path0+'nonco_angle_out.dat')
            print path0
            pprint.pprint(parser_msgs)
            #pprint.pprint(out_dict)
            cmp_values['rms'].append(out_dict['convergence_group']['rms'])
            cmp_values['charges'].append(out_dict['charge_valence_states_per_atom'])
        else:
            print path0, 'not found'
            success = False
        assert success
    pprint.pprint(cmp_values)
    # check consistency across results
    if cmp_values['rms'] != []:
        s_rms = std(cmp_values['rms'])
        pprint.pprint('std_rms= {}'.format(s_rms))
        max_s_charges = max(std(cmp_values['charges'], axis=0))
        pprint.pprint('max_std_charges= {}'.format(max_s_charges))
        assert s_rms < 10**-12
        assert max_s_charges < 10**-12

def read_file(path):
   """
   helper function to read in text file given in 'path'
   returns a flat (1D) array of numbers (floats) and list of text entries (non-cervertible to floats) of all lines that are not comments (starting with '#') 
   """
   txt = open(path).readlines()
   numbers, text = [], []
   for line in txt:
       if line[0] != '#':
          for i in line.split():
             try:
               tmp = float(i)
               numbers.append(tmp)
             except:
               text.append(i)
   return array(numbers), text

