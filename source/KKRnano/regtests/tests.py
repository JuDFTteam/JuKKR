import unittest
import subprocess as sp
import os
import glob
import shutil
import re

TESTDIR  = os.getcwd()
DECIMALS = 6

def run_it(cmd):
    """Run cmd, suppressing output.
       Returns output from stdout and exit code"""
    proc = sp.Popen(cmd, stdout=sp.PIPE, close_fds=True, preexec_fn=os.setsid, shell=True)
    out, err = proc.communicate()
    return out, proc.returncode

def get_energy(string):
    # for loop: get last match only
    match = None
    for match in re.finditer(r"^.*TOTAL ENERGY in ryd. :(.*)$", string, re.M):
        pass

    if match is not None:
         energy = float(match.group(1))
         return energy

def run_kkr(dir):
    """Run kkr-calculation with input from 'dir' and return energy"""
    setup_calc(dir)
    run_it("runkkr0.sh")
    out, err = run_it("runkkr2.sh")
    energy = get_energy(out)
    return energy

def run_kkr_mpi(dir, nranks):
    """Run kkr-calculation with input from 'dir' and return energy"""
    setup_calc(dir)
    run_it("runkkr0.sh")
    out, err = run_it("runkkr2mpi.sh {0}".format(int(nranks)))
    energy = get_energy(out)
    return energy

def setup_calc(inputdir):
    for file in glob.glob(os.path.join(inputdir, '*')):
        shutil.copy(file, TESTDIR)
 
class Test_copper(unittest.TestCase):
     def test_Cu02_lmax3(self):
        energy = run_kkr("Cu02_lmax3")
        self.assertAlmostEqual(energy, -13219.39043115, DECIMALS)

     def test_Cu02_lmax4(self):
        """Test with high lmax. Works only with -heap-arrays on ifort"""
        energy = run_kkr("Cu02_lmax4")
        self.assertAlmostEqual(energy, -13224.24988827, DECIMALS)

     def test_Cu02_lmax5(self):
        """Test with high lmax. Works only with -heap-arrays on ifort"""
        energy = run_kkr("Cu02_lmax5")
        self.assertAlmostEqual(energy, -13224.24988827, DECIMALS)

     def test_Cu02_lmax6(self):
        """Test with high lmax. Works only with -heap-arrays on ifort"""
        energy = run_kkr("Cu02_lmax6")
        self.assertAlmostEqual(energy, -13224.24988827, DECIMALS)
class Test_alloys(unittest.TestCase):

     def test_GaN(self):
        energy = run_kkr("GaN")
        self.assertAlmostEqual(energy, -13219.39043115, DECIMALS)

     def test_ZnO(self):
        energy = run_kkr("ZnO")
        self.assertAlmostEqual(energy, -13219.39043115, DECIMALS)

#     def test_Fe8Co8(self):
#        energy = run_kkr("Fe8Co8")
#        self.assertAlmostEqual(energy, -13219.39043115, DECIMALS)
class Test_Silicon(unittest.TestCase):

     def test_Si(self):
        energy = run_kkr("Si")
        self.assertAlmostEqual(energy, -13219.39043115, DECIMALS)


unittest.main()
