# coding: utf-8

"""Script to read forces from direct access file 'forces' """

import struct

DOUBLESIZE = 8
f = open('forces', 'rb')
fout = open('forces.dat', 'w')

while True:
    x = f.read(DOUBLESIZE*3)
    if x == '':
        break
    force = struct.unpack("3d", x)
    fout.write("{} {} {}\n".format(*force))

f.close()
fout.close()
