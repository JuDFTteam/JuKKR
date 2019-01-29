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
            #print "lmax = {0}".format(int(lmax))
    if solver != DEFAULT_solver:
        with open("input.conf", "a") as myfile: ## append to file
            myfile.write("solver = {0}\n".format(int(solver)))
            #print "solver = {0}".format(int(solver))
    if Lly != DEFAULT_Lly:
        with open("input.conf", "a") as myfile: ## append to file
            myfile.write("LLY = {0}\n".format(int(Lly)))
            #print "lmax = {0}".format(int(lmax))

    out, err, tim = run_it("./kkr.exe --prepare") ### start from JM-formatted potential file
    ## execute the code
    out, err, tim = run_it("OMP_STACKSIZE=80M OMP_NUM_THREADS={0} mpirun -np {1} kkr.exe".format(int(nthreads), int(nranks)))
    ### grep the result
    total_energy = get_energy(out)
    print "KKR for", inputdir, "with  lmax=",lmax,
    if solver != DEFAULT_solver:
        print ", solver=",solver,
    if nthreads != DEFAULT_nthreads:
        print ", nthreads=",nthreads,
    if nranks != DEFAULT_nranks:
        print ", nranks=",nranks,
    if Lly != DEFAULT_Lly:
        print ", Lly=",Lly,
    print " gives", total_energy, "Ryd in", tim," sec"
    out, err, tim = run_it("./clearfiles.sh")
    return total_energy

class Test_nocosocmaterials(unittest.TestCase):
    def test_MnGeB20(self):
        """Test chiral magnet MnGe B20 structure (8 atoms in unit cell)"""
        Etot = -26017.23757851
        self.assertAlmostEqual(KKRnano("MnGeB20", solver=direct, nranks=8), Etot, DECIMALS) # takes longer than other tests 
        self.assertAlmostEqual(KKRnano("MnGeB20", solver=iterative, nranks=4), Etot, DECIMALS) # takes longer than other tests 

unittest.main()
