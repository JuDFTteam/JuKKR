#!/usr/bin/env python

import pytest
from numpy import array, mean, sum, loadtxt, std

class Tests_scf_noSOC():
   """
   Tests for different mpi and hybrid runs of single iteration host_in_host
   test system: host_in_host_bccFe_bulk_FP_lmax2
   """
   def test_noSOC_mpi_serial_run(self):
      standard_scf_test('test_run_mpi_1')

   def test_noSOC_mpi_8(self):
      standard_scf_test('test_run_mpi_8')

   def test_noSOC_hybrid_1_1(self):
      standard_scf_test('test_run_hybrid_1_1')

   def test_noSOC_hybrid_1_8(self):
      standard_scf_test('test_run_hybrid_1_8')

   def test_noSOC_hybrid_8_1(self):
      standard_scf_test('test_run_hybrid_8_1')

   def test_noSOC_hybrid_2_4(self):
      standard_scf_test('test_run_hybrid_2_4')

   def test_noSOC_hybrid_4_2(self):
      standard_scf_test('test_run_hybrid_4_2')


class Tests_scf_SOC():
   """
   Check SOC scf run with different features:
    - parallelization mpi and hybrid
    - SRATRICK
    - storing wavefunctions
   """
   def test_SOC_mpi(self):
      standard_scf_test('test_run_tmatnew_mpi_8', cmpvals=[25.494781, 1.1302e-03, 3.181288, -12729.32115292])

   def test_SOC_hybrid(self):
      standard_scf_test('test_run_tmatnew_hybrid_2_4', cmpvals=[25.494781, 1.1302e-03, 3.181288, -12729.32115292])

   def test_SOC_nosavewf(self):
      standard_scf_test('test_run_nosavewf', cmpvals=[25.494781, 1.1302e-03, 3.181288, -12729.32115292])

   def test_SOC_nosratrick(self):
      standard_scf_test('test_run_nosratrick', cmpvals=[25.494780, 1.1313e-03, 3.181288, -12729.32114985])

   def test_Jij(self):
      check_Jijs('test_run_Jij', refpath='host_in_host_Jijs')

   def test_Jij_hybrid(self):
      check_Jijs('test_run_Jij_hybrid', refpath='host_in_host_Jijs')

   def test_Jij_savewf(self):
      check_Jijs('test_run_Jij_savewf', refpath='host_in_host_Jijs')

   def test_Jij_nosratrick(self):
      check_Jijs('test_run_Jij_nosratrick', refpath='host_in_host_Jijs', sracomp=True)


#helper functions

def standard_scf_test(path, cmpvals=[26., 3.1453E-08, 2.604665, -2541.14746731]):
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
   
   charge_ref = cmpvals[0]
   rms_ref    = cmpvals[1]
   spin_ref   = cmpvals[2]
   etot_ref   = cmpvals[3]
   
   assert mean(abs(charges-charge_ref)) < 10**-7
   assert abs(rms-rms_ref) < 10**-6
   assert mean(abs(spins-spin_ref)) < 10**-7
   #assert sum(abs(etot/len(charges)-etot_ref)) < 10**-6

def check_Jijs(path, refpath, sracomp=False):
   d = loadtxt('test_case_kkrflex_host_in_host/imp/'+path+'/out_Jijmatrix')
   d0= loadtxt('test_case_kkrflex_host_in_host/imp/'+refpath+'/out_Jijmatrix')

   if not sracomp:
      assert mean(abs(d-d0)) < 1e-10
      assert std(abs(d-d0)) < 1e-10
   else:
      assert mean(abs(d-d0)) < 5e-8
      assert std(abs(d-d0)) < 5e-8
