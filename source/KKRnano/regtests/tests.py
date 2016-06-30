#!/usr/bin/python

import unittest
import subprocess
import os
import glob
import shutil
import re

TESTDIR = os.getcwd() ### perform the calculation in the current working directory
DECIMALS = 6 ### 8=all digits, 6 should be enough
DEFAULT_lmax = 3
DEFAULT_solver = 3 ## solver=3 is iterative while solver=4 is direct
DEFAULT_lly = 0 ## lly=0/1 deactivate/activate Lloyd's formula

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

def KKR_total_energy(inputdir, nranks=1, nthreads=1, solver=DEFAULT_solver, lmax=DEFAULT_lmax, lly=DEFAULT_lly):
    """Run KKR-calculation with input from 'inputdir' and returns the total energy"""
    #print "start KKR for", inputdir, "with  lmax=",lmax, ", solver=",solver, ", nthreads=",nthreads, "nranks=",nranks
    out, err = run_it("./clearfiles.sh")
    
    for file in glob.glob(os.path.join(inputdir, '*')):
        shutil.copy(file, TESTDIR) ### copy all files from the input directory
        
    if lmax != DEFAULT_lmax:
        with open("input.conf", "a") as myfile: ## append to file
            myfile.write("LMAXD = {0}\n".format(int(lmax)))
            #print "lmax = {0}".format(int(lmax))
    if solver != DEFAULT_solver:
        with open("input.conf", "a") as myfile: ## append to file
            myfile.write("solver = {0}\n".format(int(solver)))
            #print "solver = {0}".format(int(solver))
    if lly != DEFAULT_lly:
        with open("input.conf", "a") as myfile: ## append to file
            myfile.write("LLY = {0}\n".format(int(lly)))
            #print "lmax = {0}".format(int(lmax))

    out, err = run_it("./kkr.exe --prepare") ### start from JM-formatted potential file
    ## execute the code
    out, err = run_it("OMP_STACKSIZE=20M OMP_NUM_THREADS={0} mpirun -np {1} kkr.exe".format(int(nthreads), int(nranks)))
    ### grep the result
    total_energy = get_energy(out)
    print "KKR for", inputdir, "with  lmax=",lmax, ", solver=",solver, ", nthreads=",nthreads, ", nranks=",nranks, " gives", total_energy, "Ryd"
    return total_energy

class Test_alloys(unittest.TestCase):
    def test_Fe8Co8(self):
       self.assertAlmostEqual(KKR_total_energy("Fe8Co8", solver=4), -42561.32515698, DECIMALS)

class Test_copper(unittest.TestCase):
     def test_Cu4_lmax(self):
        """Test with high lmax. Works only with -heap-arrays on ifort"""
        self.assertAlmostEqual(KKR_total_energy("Cu4", solver=4, lmax=3), -13219.36206406, DECIMALS)
        self.assertAlmostEqual(KKR_total_energy("Cu4", solver=4, lmax=4), -13219.71616303, DECIMALS)
        self.assertAlmostEqual(KKR_total_energy("Cu4", solver=4, lmax=5), -13219.60162033, DECIMALS)
        self.assertAlmostEqual(KKR_total_energy("Cu4", solver=4, lmax=6), -13219.56030377, DECIMALS)

class Test_semiconductors(unittest.TestCase):
     def test_GaN(self):
        self.assertAlmostEqual(KKR_total_energy("GaN", solver=4), -3990.85150060, DECIMALS)
        
     def test_Si(self):
        self.assertAlmostEqual(KKR_total_energy("Si", solver=4), -1155.68952256, DECIMALS)
        
     def test_ZnO(self):
        Etot = -7405.77074357
        self.assertAlmostEqual(KKR_total_energy("ZnO"),           Etot, DECIMALS)
        self.assertAlmostEqual(KKR_total_energy("ZnO", nranks=2), Etot, DECIMALS)
        self.assertAlmostEqual(KKR_total_energy("ZnO", nranks=4), Etot, DECIMALS)
        self.assertAlmostEqual(KKR_total_energy("ZnO", nranks=8), Etot, DECIMALS)

        Etot = -7405.77074351
        self.assertAlmostEqual(KKR_total_energy("ZnO", solver=4),           Etot, DECIMALS)
        self.assertAlmostEqual(KKR_total_energy("ZnO", solver=4, nranks=2), Etot, DECIMALS)
        self.assertAlmostEqual(KKR_total_energy("ZnO", solver=4, nranks=4), Etot, DECIMALS)
        self.assertAlmostEqual(KKR_total_energy("ZnO", solver=4, nranks=8), Etot, DECIMALS)

        Etot = -7405.74826397
        self.assertAlmostEqual(KKR_total_energy("ZnO", nranks=8, lly=1), Etot, DECIMALS)
	Etot = -7405.74826372
        self.assertAlmostEqual(KKR_total_energy("ZnO", solver=4, nranks=8, lly=1), Etot, DECIMALS)

unittest.main()
out, err = run_it("./clearfiles.sh")
