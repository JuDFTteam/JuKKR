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
        standard_verify(path0)

    def test_verify2_Fe_slab(self):
        path0 = 'test_run2_serial_1_1/'
        standard_verify(path0)

    def test_verify3_Si_lloyd(self):
        path0 = 'test_run3_serial_1_1/'
        standard_verify(path0)

    def test_compare_modes1(self):
        cmplist = ['serial_1_1', 'omp_1_1', 'omp_4_1', 'mpi_1_1' ,'mpi_1_4', 'hybrid_1_1', 'hybrid_1_4', 'hybrid_4_1', 'hybrid_2_2']
        path00 = 'test_run1_'
        cmp_values = cmp_modes(cmplist, path00)
        if cmp_values['rms'] != []:
            assert std(cmp_values['rms'])<10**-8

    def test_compare_modes2(self):
        cmplist = ['serial_1_1', 'omp_1_1', 'omp_4_1', 'mpi_1_1' ,'mpi_1_4', 'hybrid_1_1', 'hybrid_1_4', 'hybrid_4_1', 'hybrid_2_2']
        path00 = 'test_run2_'
        cmp_values = cmp_modes(cmplist, path00)
        if cmp_values['rms'] != []:
            assert std(cmp_values['rms'])<10**-8

    def test_compare_modes3(self):
        cmplist = ['serial_1_1', 'omp_1_1', 'omp_4_1', 'mpi_1_1' ,'mpi_1_4', 'hybrid_1_1', 'hybrid_1_4', 'hybrid_4_1', 'hybrid_2_2']
        path00 = 'test_run3_'
        cmp_values = cmp_modes(cmplist, path00)
        if cmp_values['rms'] != []:
            assert std(cmp_values['rms'])<10**-8

        
# helper functions

def standard_verify(path0, rms_threshold=10**-8, neutr_threshold=10**-6):
    success, parser_msgs, out_dict = parse_kkr_outputfile({}, path0+'out_kkr', path0+'output.0.txt', path0+'output.000.txt', path0+'out_timing.000.txt', path0+'out_potential', path0+'nonco_angle_out.dat')
    pprint.pprint(parser_msgs)
    pprint.pprint(out_dict)
    assert success
    assert out_dict['convergence_group']['rms'] < rms_threshold
    assert abs(out_dict['convergence_group']['charge_neutrality']) <= neutr_threshold

def cmp_modes(cmplist, path00):
    cmp_values = {'rms':[], 'charges':[]}
    for addpath in cmplist:
        path0 = path00+addpath
        print path0, os.listdir('.')
        if path0 in os.listdir('.'):
            path0+='/'
            success, parser_msgs, out_dict = parse_kkr_outputfile({}, path0+'out_kkr', path0+'output.0.txt', path0+'output.000.txt', path0+'out_timing.000.txt', path0+'out_potential', path0+'nonco_angle_out.dat')
            print path0
            pprint.pprint(parser_msgs)
            assert success
            cmp_values['rms'].append(out_dict['convergence_group']['rms'])
            cmp_values['charges'].append(out_dict['total_charge_per_atom'])
    pprint.pprint(cmp_values)
    return cmp_values  

