#!/usr/bin/python

import unittest
import subprocess
import os
import glob
import shutil
import re

TESTDIR = os.getcwd()
decimals = 6
DEFAULT_LMAX = 3

def run_it(cmd):
    """Run cmd, suppressing output. Returns output from stdout and exit code"""
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, close_fds=True, preexec_fn=os.setsid, shell=True)
    out, err = proc.communicate()
    return out, proc.returncode

def get_energy(string):
    try:
          match = list(re.finditer(r"^.*TOTAL ENERGY in ryd. :(.*)$", string, re.M))[-1] # get last match only
    except:
          print string
          raise ArgumentError
    if match is not None:
          return float(match.group(1))
    else:
          raise ArgumentError

def KKR_total_energy(dir, nranks=1, nthreads=1, method="direct", LMAX=DEFAULT_LMAX):
    """Run KKR-calculation with input from 'dir' and returns the total energy"""
    setup_calc(dir, LMAX)
    run_it("./kkr.exe --prepare")
    out, err = run_it("OMP_STACKSIZE=20M OMP_NUM_THREADS={0} mpirun -np {1} kkr.exe".format(int(nthreads), int(nranks)))
    total_energy = get_energy(out)
    return total_energy

def setup_calc(inputdir, LMAX):
    run_it("./clearfiles.sh")
    for file in glob.glob(os.path.join(inputdir, '*')):
        shutil.copy(file, TESTDIR)

class Test_semiconductors(unittest.TestCase):
     def test_ZnO(self):
        self.assertAlmostEqual(KKR_total_energy("ZnO"), -7405.77074351, decimals)

     def test_Si(self):
        self.assertAlmostEqual(KKR_total_energy("Si"), -1155.68952256, decimals)

     def test_GaN(self):
        self.assertAlmostEqual(KKR_total_energy("GaN"), -3990.85150060, decimals)

class Test_copper(unittest.TestCase):
     def test_Cu4_lmax(self):
        """Test with high lmax. Works only with -heap-arrays on ifort"""
        self.assertAlmostEqual(KKR_total_energy("Cu02_lmax3", LMAX=3), -13219.36206406, decimals)
        self.assertAlmostEqual(KKR_total_energy("Cu02_lmax4", LMAX=4), -13224.71616303, decimals)
        self.assertAlmostEqual(KKR_total_energy("Cu02_lmax5", LMAX=5), -13224.60162033, decimals)
        self.assertAlmostEqual(KKR_total_energy("Cu02_lmax6", LMAX=6), -13224.56030377, decimals)

class Test_alloys(unittest.TestCase):
    def test_Fe8Co8(self):
       self.assertAlmostEqual(KKR_total_energy("Fe8Co8"), -42561.32515698, decimals)



unittest.main()
