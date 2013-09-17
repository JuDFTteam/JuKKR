# coding: utf-8

"""Script to read forces from direct access file 'forces' """

import struct

DOUBLESIZE = 8

def forces_reader(f):
    """Reader to extract forces from binary file f."""
    while True:
        x = f.read(DOUBLESIZE*3)
        if x == '':
            break
        force = struct.unpack("3d", x)
        # Force components are not in correct order for cartesian output
        # See definition of real spherical harmonics with l=1
        yield force[2], force[0], force[1]


def main():
    """Main routine: Read forces from binary file and write to text file."""

    f = open('forces', 'rb')
    fout = open('forces.dat', 'w')
    
    force_file = forces_reader(f)
    
    for force in force_file:
        fout.write("{0} {1} {2}\n".format(*force))
    
    f.close()
    fout.close()

if __name__ == "__main__":
    main()
