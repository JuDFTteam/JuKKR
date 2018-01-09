#!/usr/bin/env python

import pytest
import os, pprint
from kkrparser_functions import parse_kkr_outputfile
from numpy import std, array


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
        cmplist = ['serial_1_1', 'omp_1_1', 'omp_4_1', 'mpi_1_1' ,'mpi_1_4', 'hybrid_1_1', 'hybrid_1_4', 'hybrid_4_1', 'hybrid_2_2']
        #cmplist = ['serial_1_1', 'omp_1_1', 'omp_4_1', 'mpi_1_1' ,'mpi_1_4', 'hybrid_1_1', 'hybrid_1_4', 'hybrid_4_1', 'hybrid_2_2']
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
        cmp_modes(cmplist, '')

        
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

