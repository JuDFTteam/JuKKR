#!/usr/bin/python

import unittest
import subprocess
import os
import glob
import shutil
import re
import time

TESTDIR = os.getcwd() ### perform the calculation in the current working directory
DECIMALS = 6 ### 8=all digits, 6 should be enough
DEFAULT_lmax = 3
DEFAULT_nranks = 1
DEFAULT_nthreads = 1
direct = 4 ##
iterative = 3 ##
DEFAULT_solver = iterative
DEFAULT_Lly = 0 ## Lly=0/1 deactivate/activate Lloyd's formula
ShowMD5 = True
AllMPIs = 1 # 1=Yes, 0=No
HighLmax = True

def run_it(cmd):
    """Run cmd, suppressing output. Returns output from stdout and exit code"""
    start_time = time.time()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, close_fds=True, preexec_fn=os.setsid, shell=True)
    out, err = proc.communicate()
    end_time = time.time()
    tim = end_time - start_time
    return out, proc.returncode, tim

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

def KKRnano(inputdir, nranks=DEFAULT_nranks, nthreads=DEFAULT_nthreads, solver=DEFAULT_solver, lmax=DEFAULT_lmax, Lly=DEFAULT_Lly):
    """Run KKR-calculation with input from 'inputdir' and returns the total energy"""
    #print "start KKR for", inputdir, "with  lmax=",lmax, ", solver=",solver, ", nthreads=",nthreads, "nranks=",nranks
    out, err, tim = run_it("./clearfiles.sh")

    global ShowMD5
    if ShowMD5:
        out, err, tim = run_it("md5sum ./kkr.exe")
        print out
        ShowMD5 = False ## do not show again
    
    for file in glob.glob(os.path.join(inputdir, '*')):
        shutil.copy(file, TESTDIR) ### copy all files from the input directory
        
    if lmax != DEFAULT_lmax:
        with open("input.conf", "a") as myfile: ## append to file
            myfile.write("LMAXD = {0}\n".format(int(lmax)))
    if solver != DEFAULT_solver:
        with open("input.conf", "a") as myfile: ## append to file
            myfile.write("solver = {0}\n".format(int(solver)))
    if Lly != DEFAULT_Lly:
        with open("input.conf", "a") as myfile: ## append to file
            myfile.write("LLY = {0}\n".format(int(Lly)))

    out, err, tim = run_it("./kkr.exe --prepare") ### start from JM-formatted potential file
    ## execute the code
    out, err, tim = run_it("OMP_STACKSIZE=80M OMP_NUM_THREADS={0} mpiexec -np {1} kkr.exe".format(int(nthreads), int(nranks)))
    ### grep the result
    total_energy = get_energy(out)
    print "KKR for",inputdir," with lmax=",lmax," gives",total_energy,"Ryd",
    if solver != DEFAULT_solver:
        print ", solver=",solver,
    if nthreads != DEFAULT_nthreads:
        print ", nthreads=",nthreads,
    if nranks != DEFAULT_nranks:
        print ", nranks=",nranks,
    if Lly != DEFAULT_Lly:
        print ", Lly=",Lly,
    print " in", tim," sec"
    out, err, tim = run_it("./clearfiles.sh")
    return total_energy

class Test_alloys(unittest.TestCase):
    def test_Fe8Co8(self):
        """Test random alloy of 16 atoms"""
        Etot = -42561.32515698
        for r in range(0, AllMPIs*4+1): # nranks=[1, 2, 4, 8, 16]
            self.assertAlmostEqual(KKRnano("Fe8Co8", solver=direct, nranks=2**r), Etot, DECIMALS) # about 30 seconds for nranks=1

class Test_copper(unittest.TestCase):
    def test_Cu4_lmax(self):
        """Test with high lmax. Works only with -heap-arrays on ifort, 4 Cu atoms in the cubic unit cell"""
        Etot = -13219.36206406
        for r in range(0, AllMPIs*2+1): # nranks=[1, 2, 4]
            self.assertAlmostEqual(KKRnano("Cu4", solver=direct, nranks=2**r), Etot, DECIMALS)
        if HighLmax:
            self.assertAlmostEqual(KKRnano("Cu4", solver=direct, lmax=4), -13219.71616303, DECIMALS)
            self.assertAlmostEqual(KKRnano("Cu4", solver=direct, lmax=5), -13219.60162033, DECIMALS) # about 30 seconds
            self.assertAlmostEqual(KKRnano("Cu4", solver=direct, lmax=6), -13219.56030377, DECIMALS) # about 60 seconds

    def test_Cu1_lmax(self):
        """Test with high lmax. Works only with -heap-arrays on ifort, 1 Cu atoms in the FCC unit cell"""
        self.assertAlmostEqual(KKRnano("Cu1", solver=direct), -3308.14107181, DECIMALS) # about  2 seconds
        if HighLmax:
            self.assertAlmostEqual(KKRnano("Cu1", solver=direct, lmax=4), -3308.26072261, DECIMALS) # about  4 seconds
            self.assertAlmostEqual(KKRnano("Cu1", solver=direct, lmax=5), -3308.22046659, DECIMALS) # about  8 seconds
            self.assertAlmostEqual(KKRnano("Cu1", solver=direct, lmax=6), -3308.15010032, DECIMALS) # about 16 seconds

class Test_semiconductors(unittest.TestCase):
    def test_GaN(self):
        """Test semiconductor in zincblende structure with 2 vacancy cells"""
        Etot = -3990.85150060
        for r in range(0, AllMPIs*2+1): # nranks=[4, 2, 1]
            self.assertAlmostEqual(KKRnano("GaN", solver=direct, nranks=2**(2-r)), Etot, DECIMALS) # about 80 seconds

    def test_Si(self):
        """Test semiconductor in diamond structure with 2 vacancy cells"""
        Etot = -1155.68952256
        for r in range(0, AllMPIs*2+1): # nranks=[1, 2, 4]
            self.assertAlmostEqual(KKRnano("Si", solver=direct, nranks=2**r), Etot, DECIMALS) # about a minute

    def test_ZnO(self):
        """Test semiconductor in wurzite structure with 4 vacancy cells and voro_weights"""
        Etot = -7405.77074357 ## test iterative solver (solver=3, default) without and with MPI
        for r in range(0, AllMPIs*3+1): # nranks=[1, 2, 4, 8]
            self.assertAlmostEqual(KKRnano("ZnO", nranks=2**r), Etot, DECIMALS)

        Etot = -7405.77074351
        for r in range(0, AllMPIs*3+1): # nranks=[1, 2, 4, 8]
            self.assertAlmostEqual(KKRnano("ZnO", solver=direct, nranks=2**r), Etot, DECIMALS)
        ### Lloyd formula
        self.assertAlmostEqual(KKRnano("ZnO", solver=direct, nranks=8, Lly=1), -7405.74826372, DECIMALS)
        self.assertAlmostEqual(KKRnano("ZnO",                nranks=8, Lly=1), -7405.74826372, DECIMALS)

unittest.main()
