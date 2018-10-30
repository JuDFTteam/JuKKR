#!/usr/bin/env python

from subprocess import check_output
from numpy import *
from matplotlib.pyplot import *


d = check_output('python test.py',shell=True).split('\n')[:-1]

data = []
for i in d:
    data.append([int(i.split()[0]),i.split()[-1]])
data = array(data)

plot(data[:,0])

d = data[:,0]
d = array([int(i) for i in d])

m = where(d[:]>10**6)

print data[m]

show()
