#!/usr/bin/env python

#use print('message') instead of print 'message' in python 2.7 as well:
from __future__ import print_function

#import pytest
import os, pprint
from kkrparser_functions import parse_kkr_outputfile
from numpy import mean, std, array, loadtxt


class Test_qdos():
    """
    check results of different features
    """
    def test_12_qdos(self):
        path  = 'test_run12_mpi_1_6/'
        path0 = 'test_run12_mpi_1_6/ref/'
        for f in 'qdos.01.1.dat qdos.01.2.dat qdos.02.1.dat qdos.02.2.dat qdos.03.1.dat qdos.03.2.dat qdos.04.1.dat qdos.04.2.dat'.split():
           fname = f
           num, text = read_file(path+fname)
           num_ref, text_ref = read_file(path0+fname)
           # remove line with serial number
           text = text[1:]
           text_ref = text_ref[1:]
           # now compare
           print(fname)
           print(std(abs(num-num_ref)))
           print(mean(abs(num-num_ref)))
           print(abs(num-num_ref).max())
           print(set(text)-set(text_ref)==set())
           assert std(abs(num-num_ref))<5*10**-16
           assert mean(abs(num-num_ref))<10**-14
           assert abs(num-num_ref).max()<2*10**-12
           assert set(text)-set(text_ref)==set()

    def test_12_qdos_SOC(self):
        path  = 'test_run12.1_mpi_1_6/'
        path0 = 'test_run12.1_mpi_1_6/ref/'
        for f in 'qdos.01.1.dat qdos.01.2.dat qdos.02.1.dat qdos.02.2.dat qdos.03.1.dat qdos.03.2.dat qdos.04.1.dat qdos.04.2.dat'.split():
           fname = f
           num, text = read_file(path+fname)
           num_ref, text_ref = read_file(path0+fname)
           # remove line with serial number
           text = text[1:]
           text_ref = text_ref[1:]
           # now compare
           print(fname)
           print(std(abs(num-num_ref)))
           print(mean(abs(num-num_ref)))
           print(abs(num-num_ref).max())
           print(set(text)-set(text_ref)==set())
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

def cmp_modes(cmplist, path00, s_rms_bound=10**-12, max_s_charges_bound=10**-12):
    """
    check convergence and charges across parallel runs
    returns dict with rms and charges entries for all paths given in 'cmplist' in parent directory 'path00'
    """
    cmp_values = {'rms':[], 'charges':[]}
    for addpath in cmplist:
        path0 = path00+addpath
        print(path0, os.listdir('.'))
        if path0 in os.listdir('.') or '/ref' in path0:
            path0+='/'
            success, parser_msgs, out_dict = parse_kkr_outputfile({}, path0+'out_kkr', path0+'output.0.txt', path0+'output.000.txt', path0+'out_timing.000.txt', path0+'out_potential', path0+'nonco_angle_out.dat')
            print(path0)
            pprint.pprint(parser_msgs)
            #pprint.pprint(out_dict)
            cmp_values['rms'].append(out_dict['convergence_group']['rms'])
            cmp_values['charges'].append(out_dict['charge_valence_states_per_atom'])
        else:
            print(path0, 'not found')
            success = False
        assert success
    pprint.pprint(cmp_values)
    # check consistency across results
    if cmp_values['rms'] != []:
        s_rms = std(cmp_values['rms'])
        pprint.pprint('std_rms= {}'.format(s_rms))
        max_s_charges = max(std(cmp_values['charges'], axis=0))
        pprint.pprint('max_std_charges= {}'.format(max_s_charges))
        assert s_rms < s_rms_bound
        assert max_s_charges < max_s_charges_bound

def read_file(path):
   """
   helper function to read in text file given in 'path'
   returns a flat (1D) array of numbers (floats) and list of text entries (non-cervertible to floats) of all lines that are not comments (starting with '#') 
   """
   txt = open(path).readlines()
   numbers, text = [], []
   for line in txt:
       # replace D+XX by e+XX for be able to read numbers in python
       line = line.replace('D', 'e').replace('\n','')
       if len(line)>0:
          if line[0] != '#':
             for i in line.split():
                try:
                  tmp = float(i)
                  numbers.append(tmp)
                except:
                  text.append(i)
   return array(numbers), text

