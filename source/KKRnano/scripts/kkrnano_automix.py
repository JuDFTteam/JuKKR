#!/usr/bin/python
import subprocess
import string
import sys


global lastline 
 
def replaceLine(file,pattern,replacement):
	text = open(file).read()
	new_text = '\n'.join(replacement if line.startswith(pattern) else line for line in text.splitlines())
	open(file, 'w').write(new_text)
 
 
def printInputParameters(mpiranks,straight_steps,anderson_steps,cycles,mixing_factor):
	print ("Number of MPI ranks: " + str(mpiranks))
	print ("Number of straight steps: " + str(straight_steps))
	print ("Number of Anderson steps: " + str(anderson_steps))
	print ("Number of Cycles: " + str(cycles))
	print ("Mixing factor: " + str(mixing_factor))
	sys.stdout.flush()
 
def mixingpattern(mpiranks,straight_steps,anderson_steps,cycles,mixing_factor):
        prepare_call = " ./kkr.exe --prepare "
	run_call = " srun -n "+str(mpiranks)+" kkr.exe "
	i = 1
	while(i<=cycles):
		if(straight_steps >= 1):
	        	print ("CONVERGING WITH STRAIGHT MIXING")
	        	replaceLine("input.conf","imix","imix = 1")
	        	replaceLine("input.conf","scfsteps","scfsteps = "+str(straight_steps))
	        	replaceLine("input.conf","mixing","mixing = "+str(mixing_factor))
			print ("Running kkr.exe --prepare...")
			sys.stdout.flush()
			subprocess.call([prepare_call], shell=True)
			print ("Restoring files from directory 'temp'... (Some error messages are expected here because not all files are needed.)")
			sys.stdout.flush()
			subprocess.call(["~/KKRnano/scripts/restore.sh temp"], shell=True)
			subprocess.call(["~/KKRnano/scripts/prepare.sh"], shell=True)
	        	print ("Running kkr.exe with straight mixing...")
			sys.stdout.flush()		
			subprocess.call([run_call], shell=True)
			subprocess.call(["~/KKRnano/scripts/save.sh temp"], shell=True)
		if(anderson_steps >= 1):
	        	print ("CONVERGING WITH ANDERSON MIXING")
	        	replaceLine("input.conf","imix","imix = 5")
	        	replaceLine("input.conf","scfsteps","scfsteps = "+str(anderson_steps))
	        	replaceLine("input.conf","mixing","mixing = "+str(mixing_factor))
			print ("Running kkr.exe --prepare...")
			sys.stdout.flush()
			subprocess.call([prepare_call], shell=True)
			print ("Restoring files from directory 'temp'... (Some error messages are expected here because not all files are needed.)")
			sys.stdout.flush()
			subprocess.call(["~/KKRnano/scripts/restore.sh temp"], shell=True)
			subprocess.call(["~/KKRnano/scripts/prepare.sh"], shell=True)
	        	print ("Running kkr.exe with Anderson mixing...")
			sys.stdout.flush()		
			subprocess.call([run_call], shell=True)
			subprocess.call(["~/KKRnano/scripts/save.sh temp"], shell=True)
		i = i + 1
     
 
# Main part, this is where the action happens
def main():
        mpiranks = int(sys.argv[1])
        straight_steps = int(sys.argv[2])
        anderson_steps = int(sys.argv[3])
        cycles = int(sys.argv[4])
        mixing_factor = float(sys.argv[5])
	printInputParameters(mpiranks,straight_steps,anderson_steps,cycles,mixing_factor)
	mixingpattern(mpiranks,straight_steps,anderson_steps,cycles,mixing_factor)
if __name__ == "__main__":
    main()
 
