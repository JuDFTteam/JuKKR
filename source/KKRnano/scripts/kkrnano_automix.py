#!/usr/bin/python
import subprocess
import string
import sys


global lastline 
 
def replaceLine(file,pattern,replacement):
	text = open(file).read()
	new_text = '\n'.join(replacement if line.startswith(pattern) else line for line in text.splitlines())
	open(file, 'w').write(new_text)
 
#def clearStringOfCharacters(text):
#	all=string.maketrans('','')
#	nodigs=all.translate(all, string.digits+"-"+".")
#	text=text.translate(all, nodigs)
#	text=text[1:]
#	return text
 
#def printTotalEnergy(file_OUT2,file_energies,alat):
#	lastline = 'NULL'
#        for line in open(file_OUT2):
#       		if 'TOTAL ENERGY' in line:			
#			lastline=line
#        print lastline
#        lastline=clearStringOfCharacters(lastline)    
#        fo = open(file_energies, "a")
#        fo.write(str((alat*0.52918)**3)+" " )
#        fo.write(str(float(lastline)*13.605698066)+"\n");
#        fo.close() 
 
#def printRMSError(file_OUT2):
#        lastline = 'NULL'
#        for line in open(file_OUT2):
#                if 'v+' in line:
#                        lastline=line
#        print lastline
 
#def printNeutrality(file_OUT2):
#        lastline = 'NULL'
#        for line in open(file_OUT2):
#                if 'neutrality' in line:
#                        lastline=line
#        print lastline
 
def printInputParameters(nodes,mpiranks,cycles,mixing_factor):
	print ("Nodes: " + str(nodes))
	print ("# MPI ranks: " + str(mpiranks))
	print ("# cycles: " + str(cycles))
	print ("Mixing factor: " + str(mixing_factor))
	sys.stdout.flush()
 
def mixingpattern(nodes,mpiranks,cycles,mixing_factor):
        prepare_call = " ./kkr.exe --prepare "
	run_call = " srun -n "+str(mpiranks)+" kkr.exe "
	i = 1
	while(i<=cycles):
		# 10 steps with straight mixing	
	        print ("1. 10 STEPS WITH STRAIGHT MIXING")
	        replaceLine("input.conf","imix","imix = 1")
	        replaceLine("input.conf","scfsteps","scfsteps = "+str(10))
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
		# 20 steps with Anderson mixing	
	        print ("1. 20 STEPS WITH Anderson MIXING")
	        replaceLine("input.conf","imix","imix = 5")
	        replaceLine("input.conf","scfsteps","scfsteps = "+str(20))
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
	nodes = int(sys.argv[1])
        mpiranks = int(sys.argv[2])
        cycles = int(sys.argv[3])
        mixing_factor = float(sys.argv[4])
	printInputParameters(nodes,mpiranks,cycles,mixing_factor)
	mixingpattern(nodes,mpiranks,cycles,mixing_factor)
if __name__ == "__main__":
    main()
 
