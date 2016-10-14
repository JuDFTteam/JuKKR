#!/usr/bin/env python

from numpy import *
from matplotlib.pyplot import *

d = loadtxt('out')

figure()
for i in range(1,4):
   m  = where(d[:,1]==i)
   tmp = d[m]
   #subplot(3,1,i)
   plot(tmp[:,-2],tmp[:,-1],'o:')

show()
