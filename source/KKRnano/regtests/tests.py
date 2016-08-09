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
DEFAULT_solver = 3 ## solver=3 is iterative while solver=4 is direct
DEFAULT_lly = 0 ## lly=0/1 deactivate/activate Lloyd's formula
ShowMD5 = True
AllMPIs = True
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

def KKRnano(inputdir, nranks=1, nthreads=1, solver=DEFAULT_solver, lmax=DEFAULT_lmax, lly=DEFAULT_lly):
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
            #print "lmax = {0}".format(int(lmax))
    if solver != DEFAULT_solver:
        with open("input.conf", "a") as myfile: ## append to file
            myfile.write("solver = {0}\n".format(int(solver)))
            #print "solver = {0}".format(int(solver))
    if lly != DEFAULT_lly:
        with open("input.conf", "a") as myfile: ## append to file
            myfile.write("LLY = {0}\n".format(int(lly)))
            #print "lmax = {0}".format(int(lmax))

    out, err, tim = run_it("./kkr.exe --prepare") ### start from JM-formatted potential file
    ## execute the code
    out, err, tim = run_it("OMP_STACKSIZE=20M OMP_NUM_THREADS={0} mpirun -np {1} kkr.exe".format(int(nthreads), int(nranks)))
    ### grep the result
    total_energy = get_energy(out)
    print "KKR for", inputdir, "with  lmax=",lmax, ", solver=",solver, ", nthreads=",nthreads, ", nranks=",nranks, " gives", total_energy, "Ryd in", tim," sec"
    out, err, tim = run_it("./clearfiles.sh")
    return total_energy

class Test_alloys(unittest.TestCase):
    def test_Fe8Co8(self):
       """Test random alloy of 16 atoms"""
       Etot = -42561.32515698
       self.assertAlmostEqual(KKRnano("Fe8Co8", solver=4), Etot, DECIMALS) # about 30 seconds
       if AllMPIs:
           self.assertAlmostEqual(KKRnano("Fe8Co8", solver=4, nranks=2), Etot, DECIMALS)
           self.assertAlmostEqual(KKRnano("Fe8Co8", solver=4, nranks=4), Etot, DECIMALS)
           self.assertAlmostEqual(KKRnano("Fe8Co8", solver=4, nranks=8), Etot, DECIMALS)
           self.assertAlmostEqual(KKRnano("Fe8Co8", solver=4, nranks=16),Etot, DECIMALS)

class Test_copper(unittest.TestCase):
    def test_Cu4_lmax(self):
        """Test with high lmax. Works only with -heap-arrays on ifort, 4 Cu atoms in the cubic unit cell"""
        Etot = -13219.36206406
        self.assertAlmostEqual(KKRnano("Cu4", solver=4), Etot, DECIMALS)
        if AllMPIs:
            self.assertAlmostEqual(KKRnano("Cu4", solver=4, nranks=2), Etot, DECIMALS)
            self.assertAlmostEqual(KKRnano("Cu4", solver=4, nranks=4), Etot, DECIMALS)
        if HighLmax:
            self.assertAlmostEqual(KKRnano("Cu4", solver=4, lmax=4), -13219.71616303, DECIMALS)
            self.assertAlmostEqual(KKRnano("Cu4", solver=4, lmax=5), -13219.60162033, DECIMALS) # about 30 seconds
            self.assertAlmostEqual(KKRnano("Cu4", solver=4, lmax=6), -13219.56030377, DECIMALS) # about 60 seconds

    def test_Cu1_lmax(self):
        """Test with high lmax. Works only with -heap-arrays on ifort, 1 Cu atoms in the FCC unit cell"""
        self.assertAlmostEqual(KKRnano("Cu1", solver=4), -3308.14107181, DECIMALS) # about  2 seconds
        if HighLmax:
            self.assertAlmostEqual(KKRnano("Cu1", solver=4, lmax=4), -3308.26072261, DECIMALS) # about  4 seconds
            self.assertAlmostEqual(KKRnano("Cu1", solver=4, lmax=5), -3308.22046659, DECIMALS) # about  8 seconds
            self.assertAlmostEqual(KKRnano("Cu1", solver=4, lmax=6), -3308.15010032, DECIMALS) # about 16 seconds

class Test_semiconductors(unittest.TestCase):
    def test_GaN(self):
        """Test semiconductor in zincblende structure with 2 vacancy cells"""
        Etot = -3990.85150060
        self.assertAlmostEqual(KKRnano("GaN", solver=4, nranks=4), Etot, DECIMALS) # about 80 seconds
        if AllMPIs:
            self.assertAlmostEqual(KKRnano("GaN", solver=4, nranks=2), Etot, DECIMALS)
            self.assertAlmostEqual(KKRnano("GaN", solver=4, nranks=1), Etot, DECIMALS) # about 4 minutes

    def test_Si(self):
        """Test semiconductor in diamond structure with 2 vacancy cells"""
        Etot = -1155.68952256
        self.assertAlmostEqual(KKRnano("Si", solver=4), Etot, DECIMALS) # about a minute
        if AllMPIs:
            self.assertAlmostEqual(KKRnano("Si", solver=4, nranks=2), Etot, DECIMALS)
            self.assertAlmostEqual(KKRnano("Si", solver=4, nranks=4), Etot, DECIMALS)

    def test_ZnO(self):
        """Test semiconductor in wurzite structure with 4 vacancy cells and voro_weights"""
        Etot = -7405.77074357 ## test iterative solver (solver=3, default) without and with MPI
        self.assertAlmostEqual(KKRnano("ZnO"), Etot, DECIMALS)
        if AllMPIs:
            self.assertAlmostEqual(KKRnano("ZnO", nranks=2), Etot, DECIMALS)
            self.assertAlmostEqual(KKRnano("ZnO", nranks=4), Etot, DECIMALS)
            self.assertAlmostEqual(KKRnano("ZnO", nranks=8), Etot, DECIMALS)

        Etot = -7405.77074351
        self.assertAlmostEqual(KKRnano("ZnO", solver=4), Etot, DECIMALS)
        if AllMPIs:
            self.assertAlmostEqual(KKRnano("ZnO", solver=4, nranks=2), Etot, DECIMALS)
            self.assertAlmostEqual(KKRnano("ZnO", solver=4, nranks=4), Etot, DECIMALS)
            self.assertAlmostEqual(KKRnano("ZnO", solver=4, nranks=8), Etot, DECIMALS)
        # Lloyd formula
        self.assertAlmostEqual(KKRnano("ZnO", solver=4, nranks=8, lly=1), -7405.74826372, DECIMALS)
        self.assertAlmostEqual(KKRnano("ZnO", nranks=8, lly=1), -7405.74826397, DECIMALS)

unittest.main()
