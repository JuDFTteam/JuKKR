#!/usr/bin/python

"""
Repeat a unit cell in each direction by given number of times
author: Elias Rabel
"""

import os
import sys

class Bravais:
    def __init__(self):
        self.bravais = []
    def getFromInputcard(self, filename):
        """Read bravais vectors from inputcard"""
        with open(filename) as file:
            for line in file:
                if line.strip() == 'BRAVAIS':
                    break   

            for ind in range(3):
                line = file.next()
                splitted = line.split()
                coords = [float(splitted[0]), float(splitted[1]), float(splitted[2])]
                self.bravais.append(coords)

    def repeat(self, times):
        rep_bravais = Bravais()
        for vec, n in zip(self.bravais, times):
            rep_bravais.bravais.append( [c*n for c in vec] )
        return rep_bravais             

class BasisCoords:
    def __init__(self):
        self.coords  = []

    def read(self, filename):
        """Read rbasis file"""
        with open(filename) as file:
	    for i in xrange(2):
	    	file.next()
            for line in file:
                splitted = line.split()
                coords = [str(splitted[0]),float(splitted[1]), float(splitted[2]), float(splitted[3])]
                self.coords.append(coords)

    def write(self, filename, num_atoms):
	     """Write rbasis file"""
	     with open(filename, "w") as file:
		     file.write(str(num_atoms) + '\n' )
		     file.write('#comment line' + '\n' )
		     for coord in self.coords:
			     file.write('     '.join( [str(c) for c in coord] ) + '\n' )

    def repeat(self, bravais, times):
        repeatedRbasis = BasisCoords()
        vectors = bravais.bravais
        repeatedRbasis.coords = []
        for iz in range(times[2]):
            for iy in range(times[1]):
               for ix in range(times[0]):
                   for coord in self.coords:
                       newcoord = []
		       newcoord.append(coord[0])
                       for ind in range(1,4):
                           newcoord.append(coord[ind] + ix * vectors[0][ind-1]
                                                      + iy * vectors[1][ind-1]
                                                      + iz * vectors[2][ind-1])
                       repeatedRbasis.coords.append(newcoord)

        return repeatedRbasis

def repeatTextFile(infilename, outfilename, number):
    """Concatenate a textfile 'number' of times"""
    with open(outfilename, "w") as outfile:
        for ind in range(number):
            with open(infilename, "r") as infile:
                for line in infile:
                    outfile.write(line)


if __name__ == "__main__":
    PATH = 'newcell'
    times = [1, 1, 1]

    num_atoms = int( sys.argv[1] )
    for ind in range(1,4):
        times[ind-1] = int( sys.argv[ind+1] )
        
    number = times[0]*times[1]*times[2]
    print("Repeat the cell " + str(number) + " times.") 
           
    rbasis = BasisCoords()
    rbasis.read('rbasis.xyz')
    print rbasis.coords

    brav = Bravais()
    brav.getFromInputcard('inputcard')
    print brav.bravais

    repbasis = rbasis.repeat(brav, times)

    print("New Bravais vectors:")
    print(brav.repeat(times).bravais)

    try:
        os.mkdir(PATH)
    except OSError:
        print("Directory " + PATH + " exists.")
        
    repbasis.write(os.path.join(PATH, 'rbasis.xyz'), num_atoms)
#    repeatTextFile('atominfo', os.path.join(PATH, 'atominfo'), number)
    repeatTextFile('potential', os.path.join(PATH, 'potential'), number)
    if os.path.isfile('nonco_angle_out.dat'):
	    repeatTextFile('nonco_angle_out.dat', os.path.join(PATH, 'nonco_angle.dat'), number)
    if os.path.isfile('voro_weights'):
	    repeatTextFile('voro_weights', os.path.join(PATH, 'voro_weights'), number)
    
