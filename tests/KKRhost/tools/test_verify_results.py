#!/usr/bin/env python

#use print('message') instead of print 'message' in python 2.7 as well:
from __future__ import print_function

#import pytest
import os, pprint
from kkrparser_functions import parse_kkr_outputfile
from numpy import mean, std, array, loadtxt


class Test_serial():
    """
    check results of simple scf tests without SOC
    """
    def test_1_Au_bulk(self):
        path0 = 'test_run01_serial_1_1/'
        standard_verify(path0, rms_threshold=5*10**-8)

    def test_2_Fe_slab(self):
        path0 = 'test_run02_serial_1_1/'
        standard_verify(path0, rms_threshold=7*10**-8, debug=True)

    def test_3_Si_lloyd(self):
        path0 = 'test_run03_serial_1_1/'
        standard_verify(path0, rms_threshold=8*10**-9)


class Test_features():
    """
    check results of different features
    """
    def test_4_Jijs_noSOC(self):
        paths = 'test_run04_serial_1_1/ test_run04_hybrid_1_3/'.split()
        path0 = 'test_inputs/test_04_Jijs_Fe_slab_lmax2_noSOC/ref/'
        # compare Jij.atom* files of different runs with reference (in path0)
        for path in paths:
           files = 'Jij.atom00002 Jij.atom00003'.split()
           for fname in files:
              num, text = read_file(path+fname)
              num_ref, text_ref = read_file(path0+fname)
              assert std(num-num_ref)<10**-10
              assert set(text)-set(text_ref)==set()

    def test_5_kkrflex(self):
        paths = 'test_run05_serial_1_1/ test_run05_hybrid_1_3/'.split()
        path0 = 'test_inputs/test_05_Silicon_lloyd_kkrflex_output_lmax2_noSOC/ref/'
        # compare kkrflex_* files of different runs with reference (in path0)
        for path in paths:
           files = 'kkrflex_atominfo kkrflex_hoststructure.dat kkrflex_intercell_cmoms kkrflex_intercell_ref kkrflex_tmat'.split()
           for fname in files:
              num, text = read_file(path+fname)
              num_ref, text_ref = read_file(path0+fname)
              assert std(num-num_ref)<10**-10
              assert set(text)-set(text_ref)==set()

    def test_6_FERMIOUT(self):
        paths = 'test_run06_serial_1_1/ test_run06_hybrid_1_3/'.split()
        path0 = 'test_inputs/test_06_Silicon_lloyd_FERMIOUT_output_lmax2_noSOC/ref/'
        # compare TBkkr_* files of different runs with reference (in path0)
        for path in paths:
           files = 'TBkkr_container.txt  TBkkr_params.txt'.split()
           for fname in files:
              num, text = read_file(path+fname)
              num_ref, text_ref = read_file(path0+fname)
              assert std(num-num_ref)<10**-10
              assert set(text)-set(text_ref)==set()

    def test_10_OPERATOR(self):
        path  = 'test_run10_mpi_1_8/'
        path0 = 'test_run10_mpi_1_8/ref/'
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
          print(filename)
          d = loadtxt(path+filename)
          d0 = loadtxt(path+filename.replace('.txt', '_imp.txt'))
          nsigma = 3
          if 'rhod' in filename:
             nsigma +=1
          # real part of middle part (should be equivalent to 'fake' impurity calculation)
          d1 = d[:,0].reshape(nsigma,22, 18, 18); d1 = d1[:,6:16,:,:]
          d01 = d0[:,0].reshape(nsigma,10,18,18)
          # imaginary part
          d2 = d[:,1].reshape(nsigma,22, 18, 18); d2 = d2[:,6:16,:,:]
          d02 = d0[:,1].reshape(nsigma,10,18,18)
          # flatten arrays and take diff
          d1 = d1.reshape(-1); d2 = d2.reshape(-1); d01 = d01.reshape(-1); d02 = d02.reshape(-1)
          diff1 = d01-d1; diff2 = d02-d2
          assert mean(diff1) < 10**-15
          assert abs(diff1).max() < 10**-15
          assert mean(diff2) < 10**-15
          assert abs(diff2).max() < 10**-15

    def test_11_DTM_GMAT(self):
        path  = 'test_run11_mpi_1_8/'
        path0 = 'test_run11_mpi_1_8/ref/'
        for f in 'DTM/DTMTRX ./green_host GMAT/GMATLL_GES'.split():
           fname = f
           num, text = read_file(path+fname)
           num_ref, text_ref = read_file(path0+fname.split('/')[1])
           print(fname)
           print(std(abs(num-num_ref)))
           print(mean(abs(num-num_ref)))
           print(abs(num-num_ref).max())
           print(set(text)-set(text_ref)==set())
           assert std(abs(num-num_ref))<5*10**-11
           assert mean(abs(num-num_ref))<10**-12
           assert abs(num-num_ref).max()<2*10**-8
           assert set(text)-set(text_ref)==set()

    """
    def test_13_rhoq(self):
        path  = 'test_run13/'
        path0 = 'test_run13/ref/'
        fname = 'out_rhoq.txt'
        num, text = read_file(path+fname)
        num_ref, text_ref = read_file(path0+fname)
        # now compare
        print(fname)
        print(std(abs(num-num_ref)))
        print(mean(abs(num-num_ref)))
        print(abs(num-num_ref).max())
        print(set(text)-set(text_ref)==set())
        assert std(abs(num-num_ref))<5*10**-9
        assert mean(abs(num-num_ref))<10**-10
        assert abs(num-num_ref).max()<2*10**-7
        assert set(text)-set(text_ref)==set()
    """

    def test_14_ASA(self):
        path0 = 'test_run14_hybrid_1_3/'
        standard_verify(path0, rms_threshold=8*10**-8, rms_threshold_end=8*10**-8, neutr_threshold=2.5*10**-5)

    def test_15_CPA(self):
        path0 = 'test_run15_hybrid_1_3/'
        standard_verify(path0, rms_threshold=8*10**-8, rms_threshold_end=8*10**-8)

    def test_16_Dirac(self):
        path0 = 'test_run16_hybrid_1_1/'
        standard_verify(path0, rms_threshold=7*10**-8, rms_threshold_end=7*10**-8)

    def test_17_lambda_xc(self):
        path0 = 'test_run17_hybrid_1_3/'
        standard_verify(path0, rms_threshold=9*10**-8, rms_threshold_end=9*10**-8)

    # no working test case yet
    #def test_18_noco(self):
    #    path0 = 'test_run18_hybrid_1_3/'
    #    standard_verify(path0, rms_threshold=9*10**-8, rms_threshold_end=9*10**-8)

    def test_19_decimate(self):
        path00 = 'test_run19_mpi_2_4'
        # first check decifile generation
        path = path00+'/bulk/'
        path0 = path+'/ref/'
        fname = 'decifile'
        num, text = read_file(path+fname)
        num_ref, text_ref = read_file(path0+fname)
        assert std(num-num_ref)<10**-10
        assert set(text)-set(text_ref)==set()
        # now check decimation step (diff against ref)
        # deactivate following test since decimation test is only valid with slab inversion (full inv gives different result)
        #cmplist = [path00, path00+'/ref']
        #cmp_modes(cmplist, '', s_rms_bound=2.13*10**-1, max_s_charges_bound=10**-3)

    def test_20_godfrin(self):
        path00 = 'test_run20_hybrid_1_3'
        # first check the two runs individually
        path = path00+'/godfrinON/'
        #standard_verify(path, rms_threshold=9*10**-8, rms_threshold_end=9*10**-8)
        # then check no godfrin run
        path0 = path00+'/godfrinOFF/'
        #standard_verify(path0, rms_threshold=9*10**-8, rms_threshold_end=9*10**-8)
        # and full inv run
        path00 = path00+'/fullinv/'
        #standard_verify(path00, rms_threshold=9*10**-8, rms_threshold_end=9*10**-8)
        # then check if output is the same for the three runs
        fname = 'output.000.txt'
        num, text = read_file(path+fname)
        # check against no godfrin
        num_ref, text_ref = read_file(path0+fname)
        assert std(num-num_ref)<10**-14
        # check against full inv
        num_ref, text_ref = read_file(path00+fname)
        assert std(num-num_ref)<2*10**-11

    def test_21_XCs(self):
        path00 = 'test_run21_hybrid_1_3'
        # first check the two runs individually
        for padd in 'MJW vBH VWN PW91 PBE PBEsol'.split():
           path = path00+'/'+padd+'/'
           standard_verify(path, rms_threshold=5*10**-7, rms_threshold_end=10**-7)

    def test_22_LDAU(self):
        path00 = 'test_run22_hybrid_1_3'
        # first check the two runs individually
        for padd in 'noSOC SOC'.split():
           path = path00+'/'+padd+'/'
           standard_verify(path, rms_threshold=5*10**-8, rms_threshold_end=5*10**-8)

    def test_23_DOS(self):
        path0 = 'test_run23_hybrid_1_3/bulk_SOC/'
        path00 = 'test_run23_hybrid_1_3/slab_noSOC/'
        for fname in 'complex.dos  dos.atom1  dos.atom2  dos.atom3  dos.atom4'.split():
           # check bulk run with SOC
           num, text = read_file(path0+fname)
           num_ref, text_ref = read_file(path0+'/ref/'+fname)
           assert std(num-num_ref)<10**-13
           # check slab run without SOC
           num, text = read_file(path00+fname)
           num_ref, text_ref = read_file(path00+'/ref/'+fname)
           assert std(num-num_ref)<10**-13


class Test_SOC():
    """
    check results of different features with SOC
    """
    def test_1_Au_bulk(self):
        path0 = 'test_run01.1_hybrid_1_3/'
        standard_verify(path0, rms_threshold=1*10**-8)

    def test_2_Fe_slab(self):
        path0 = 'test_run02.1_hybrid_1_3/'
        standard_verify(path0, rms_threshold=1*10**-8, rms_threshold_end=1*10**-8)

    def test_3_Si_lloyd(self):
        path0 = 'test_run03.1_hybrid_1_3/'
        standard_verify(path0, rms_threshold=2*10**-8, rms_threshold_end=2*10**-8)

    def test_3_2_NOSOC(self):
        path0 = 'test_run03.2_hybrid_1_3/'
        # check convergence of both runs
        standard_verify(path0+'NEWSOSOL_NOSOC/', rms_threshold=1.5*10**-6, rms_threshold_end=1.5*10**-6, neutr_threshold=8*10**-5)
        standard_verify(path0+'NEWSOSOL_SOCSCL0/', rms_threshold=1.5*10**-6, rms_threshold_end=1.5*10**-6, neutr_threshold=8*10**-5)
        # cross check both runs against each other (comparing output writte to 'out_last.txt')
        num, text = read_file(path0+'NEWSOSOL_NOSOC/out_last.txt')
        num_ref, text_ref = read_file(path0+'NEWSOSOL_SOCSCL0/out_last.txt')
        assert std(num-num_ref)<5*10**-12
        assert set(text)-set(text_ref)==set()

    def test_4_Jijs_SOC(self):
        paths = ['test_run04.1_hybrid_1_3/']
        path0 = 'test_inputs/test_04.1/ref/'
        # compare Jij.atom* files of different runs with reference (in path0)
        for path in paths:
           files = 'Jij.atom00002 Jij.atom00003'.split()
           for fname in files:
              num, text = read_file(path+fname)
              num_ref, text_ref = read_file(path0+fname)
              assert std(num-num_ref)<10**-10
              assert set(text)-set(text_ref)==set()

    def test_5_kkrflex(self):
        paths = ['test_run05.1_hybrid_1_3/']
        path0 = 'test_inputs/test_05.1/ref/'
        # compare kkrflex_* files of different runs with reference (in path0)
        for path in paths:
           files = 'kkrflex_atominfo kkrflex_hoststructure.dat kkrflex_intercell_cmoms kkrflex_intercell_ref kkrflex_tmat'.split()
           for fname in files:
              num, text = read_file(path+fname)
              num_ref, text_ref = read_file(path0+fname)
              assert std(num-num_ref)<2*10**-8
              assert set(text)-set(text_ref)==set()

    def test_6_FERMIOUT(self):
        paths = ['test_run06.1_hybrid_1_3/']
        path0 = 'test_inputs/test_06.1/ref/'
        # compare TBkkr_* files of different runs with reference (in path0)
        for path in paths:
           files = 'TBkkr_container.txt  TBkkr_params.txt'.split()
           for fname in files:
              num, text = read_file(path+fname)
              num_ref, text_ref = read_file(path0+fname)
              assert std(num-num_ref)<10**-10
              assert set(text)-set(text_ref)==set()

    def test_14_ASA(self):
        path0 = 'test_run14.1_hybrid_1_3/'
        standard_verify(path0, rms_threshold=3*10**-8, rms_threshold_end=3*10**-8, neutr_threshold=1.5*10**-5)

    def test_15_CPA(self):
        path0 = 'test_run15.1_hybrid_1_3/'
        standard_verify(path0, rms_threshold=5*10**-8, rms_threshold_end=5*10**-8)

    def test_19_decimate(self):
        path00 = 'test_run19.1_mpi_2_4'
        # first check decifile generation
        path = path00+'/bulk/'
        path0 = path+'ref/'
        fname = 'decifile'
        num, text = read_file(path+fname)
        num_ref, text_ref = read_file(path0+fname)
        assert std(num-num_ref)<10**-10
        assert set(text)-set(text_ref)==set()
        # deactivate following test since decimation test is only valid with slab inversion (full inv gives different result)
        ## now check decimation step (diff against ref)
        #cmplist = [path00, path00+'/ref']
        #cmp_modes(cmplist, '', s_rms_bound=3.28*10**-1, max_s_charges_bound=10**-3)
        
# helper functions

def standard_verify(path0, rms_threshold=10**-8, rms_threshold_end=10**-8, neutr_threshold=10**-6, debug=False):
    """
    wrapper for standard tests reading output and comparins rms and charge neutrality
    """
    # use parser function from aiida-kkr
    success, parser_msgs, out_dict = parse_kkr_outputfile({}, path0+'out_kkr', path0+'output.0.txt', path0+'output.000.txt', path0+'out_timing.000.txt', path0+'out_potential', path0+'nonco_angle_out.dat', debug=debug)
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

