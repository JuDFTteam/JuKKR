#!/usr/bin/env python

import pytest

class Tests_scf_host_in_host_bccFe_bulk_FP_noSOC_lmax2():
   """
   tests for different mpi and hybrid runs of single iteration host_in_host
   """
   def test_host_in_host_mpi_serial_run(self):
      standard_scf_test('test_run_mpi_1')

   def test_host_in_host_mpi_8(self):
      standard_scf_test('test_run_mpi_8')

   def test_host_in_host_hybrid_1_1(self):
      standard_scf_test('test_run_hybrid_1_1')

   def test_host_in_host_hybrid_1_8(self):
      standard_scf_test('test_run_hybrid_1_8')

   def test_host_in_host_hybrid_8_1(self):
      standard_scf_test('test_run_hybrid_8_1')

   def test_host_in_host_hybrid_2_4(self):
      standard_scf_test('test_run_hybrid_2_4')

   def test_host_in_host_hybrid_4_2(self):
      standard_scf_test('test_run_hybrid_4_2')
      


#helper functions

def standard_scf_test(path):
   from numpy import array, mean, sum
   d = open('test_case_kkrflex_host_in_host/imp/'+path+'/out').readlines()
   charges, spins = [], []
   for i in d:
      if 'charge in WS-cell =' in i:
          ch = float(i.split()[6])
          s = float(i.split()[10])
          charges.append(ch)
          spins.append(s)
      if 'average rms-error' in i:
          rms = float(i.split()[-1].replace('D','E'))
      if 'TOTAL ENERGY in ryd' in i:
          etot = float(i.split()[-1])
   charges = array(charges)
   spins = array(spins)
   
   charge_ref = 26.
   rms_ref = 3.1453E-08
   spin_ref = 2.604665 
   etot_ref = -2541.14746731 
   
   assert mean(abs(charges-charge_ref)) < 10**-7
   assert abs(rms-rms_ref) < 10**-6
   assert mean(abs(spins-spin_ref)) < 10**-7
   #assert sum(abs(etot/len(charges)-etot_ref)) < 10**-6
