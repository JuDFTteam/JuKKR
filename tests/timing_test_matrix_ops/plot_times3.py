#!/usr/bin/env python
from numpy import *
from matplotlib.pyplot import *
import os

if 'out' in os.listdir('.'):
  d = loadtxt('out')

  subplot(3,2,1)
  plot(d[::3,0], d[::3,2], label='slice'); plot(d[::3,0], d[::3,3], label='direct loop'); plot(d[::3,0], d[::3,1], label='zcopy')
  legend(); xlabel('matrix size'); ylabel('time for 10000 in seconds')
  subplot(3,2,2)
  plot(d[::3,0], d[::3,2]/d[::3,1], label='slice/zcopy'); plot(d[::3,0], d[::3,3]/d[::3,1], label='direct/zcopy')
  legend(loc=4); xlabel('matrix size')
  subplot(3,2,3)
  plot(d[1::3,0], d[1::3,2], label='slice'); plot(d[1::3,0], d[1::3,3], label='direct loop'); plot(d[1::3,0], d[1::3,1], label='zscal')
  legend(); xlabel('matrix size'); ylabel('time for 10000 in seconds')
  subplot(3,2,4)
  plot(d[1::3,0], d[1::3,2]/d[1::3,1], label='slice/zscal'); plot(d[1::3,0], d[1::3,3]/d[1::3,1], label='direct/zscal')
  legend(loc=4); xlabel('matrix size')
  subplot(3,2,5)
  plot(d[2::3,0], d[2::3,2], label='slice'); plot(d[2::3,0], d[2::3,3], label='direct loop'); plot(d[2::3,0], d[2::3,1], label='zaxpy')
  legend(); xlabel('matrix size'); ylabel('time for 10000 in seconds')
  subplot(3,2,6)
  plot(d[2::3,0], d[2::3,2]/d[2::3,1], label='slice/zaxpy'); plot(d[2::3,0], d[2::3,3]/d[2::3,1], label='direct/zaxpy')
  legend(loc=4); xlabel('matrix size')

  show()
