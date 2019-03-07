#!/usr/bin/env python

import sys
from numpy import loadtxt, sqrt, sum, log, abs
from matplotlib.pyplot import imshow, colorbar, show, figure, subplot, axis, subplots_adjust, title, plot, legend, gca, axvline

if len(sys.argv)>1:
  ir = int(sys.argv[1])
else:
  ir = 0
if len(sys.argv)>2:
  lin = False
else:
  lin = True

r = loadtxt('test_rmesh')
rpan = loadtxt('test_rpan_intervall')

def plot_wavefun(name, ir, lin, tadd, ipl):
  d = loadtxt('test_'+name)

  figure(4+icomp)
  subplot(2,2,1+ipl)
  title(name+tadd)
  l0 = int(sqrt(len(d)))
  d1 = d.reshape(l0,l0,-1)
  if lin:
    plot(r, d1[0, 0, 2+icomp::2], label='(lm1,lm2)=(1,1)')
    plot(r, d1[1, 1, 2+icomp::2], label='(lm1,lm2)=(2,2)')
    plot(r, d1[4, 4, 2+icomp::2], label='(lm1,lm2)=(5,5)')
  else:
    plot(r, abs(d1[0, 0, 2+icomp::2]), label='(lm1,lm2)=(1,1)')
    plot(r, abs(d1[1, 1, 2+icomp::2]), label='(lm1,lm2)=(2,2)')
    plot(r, abs(d1[4, 4, 2+icomp::2]), label='(lm1,lm2)=(5,5)')
    ax = gca()
    ax.set_yscale('log')
  if ir!=0:
    axvline(r[abs(ir)], color='r', lw=0.5)
  legend()

  figure(2+icomp)
  subplot(2,2,1+ipl)
  title(name+tadd)
  if ir<=0:
    d = sum(d[:,2+icomp+ir::2], axis=1)
  else:
    d = d[:,2+ir+icomp]
  d = d.reshape(int(sqrt(len(d))),-1)

  if lin:
    imshow(d, interpolation='nearest')
  else:
    imshow(log(abs(d)), interpolation='nearest')

  axis('equal')
  colorbar()


#figsize=(12,6))
#subplots_adjust(left=0.04, bottom=0.25, right=0.99, top=0.75)

# first kind
figure(100)
d = loadtxt('test_wronskian')
l0 = int(sqrt(len(d)))
d1 = d.reshape(l0,l0,-1)
icomp = 0
if lin:
  plot(r, d1[0, 0, 2+icomp::2], label='(lm1,lm2)=(1,1)')
  plot(r, d1[1, 1, 2+icomp::2], label='(lm1,lm2)=(2,2)')
  plot(r, d1[4, 4, 2+icomp::2], label='(lm1,lm2)=(5,5)')
else:
  #plot(r, abs(d1[0, 0, 2+icomp::2]-d1[0, 0, 2+icomp::2][-1]), label='(lm1,lm2)=(1,1)')
  #plot(r, abs(d1[1, 1, 2+icomp::2]-d1[1, 1, 2+icomp::2][-1]), label='(lm1,lm2)=(2,2)')
  #plot(r, abs(d1[4, 4, 2+icomp::2]-d1[4, 4, 2+icomp::2][-1]), label='(lm1,lm2)=(5,5)')
  plot(r, abs(d1[0, 0, 2+icomp::2]), label='(lm1,lm2)=(1,1)')
  plot(r, abs(d1[1, 1, 2+icomp::2]), label='(lm1,lm2)=(2,2)')
  plot(r, abs(d1[4, 4, 2+icomp::2]), label='(lm1,lm2)=(5,5)')
  ax = gca()
  ax.set_xlim(r[1], r[-1])
  ax.set_yscale('log')
  ax.set_xscale('log')
if ir!=0:
  axvline(r[abs(ir)], color='r', lw=0.8)
for ipan in rpan:
  axvline(ipan, color='grey', lw=0.5)
"""
ax2 = gca().twiny()
ax2.plot(range(len(r))) # Create a dummy plot
ax2.cla()
if not lin:
  ax2.set_xlim(1, len(r)-1)
  ax2.set_xscale('log')
"""

figure(1)
subplot(2,2,1)
title('wronskian of first kind')
if ir<=0:
  d = sum(d[:,2+ir::2], axis=1)
else:
  d = d[:,2+ir]
d = d.reshape(int(sqrt(len(d))),-1)
if lin:
  imshow(d, interpolation='nearest')
else:
  imshow(log(abs(d)), interpolation='nearest')
axis('equal')
colorbar()
# imag
figure(1)
subplot(2,2,2)
title('imaginary part')
if ir<=0:
  d = sum(loadtxt('test_wronskian')[:,3+ir::2], axis=1)
else:
  d = loadtxt('test_wronskian')[:,3+ir]
d = d.reshape(int(sqrt(len(d))),-1)
if lin:
  imshow(d, interpolation='nearest')
else:
  imshow(log(abs(d)), interpolation='nearest')
axis('equal')
colorbar()

# second kind
figure(1)
subplot(2,2,3)
title('wronskian of second kind')
if ir<=0:
  d = sum(loadtxt('test_wronskian2')[:,2+ir::2], axis=1)
else:
  d = loadtxt('test_wronskian2')[:,2+ir]
d = d.reshape(int(sqrt(len(d))),-1)
if lin:
  imshow(d, interpolation='nearest')
else:
  imshow(log(abs(d)), interpolation='nearest')
axis('equal')
colorbar()
# imag
figure(1)
subplot(2,2,4)
title('imaginary part')
if ir<=0:
  d = sum(loadtxt('test_wronskian2')[:,3+ir::2], axis=1)
else:
  d = loadtxt('test_wronskian2')[:,3+ir]
d = d.reshape(int(sqrt(len(d))),-1)
if lin:
  imshow(d, interpolation='nearest')
else:
  imshow(log(abs(d)), interpolation='nearest')
axis('equal')
colorbar()

#"""
for icomp in [0,1]: # for real/imag
  if icomp==0:
    tadd = ' real part'
  else:
    tadd = ' imag part'
  plot_wavefun('rll', ir, lin, tadd, 0)
  plot_wavefun('rllleft', ir, lin, tadd, 1)
  plot_wavefun('sll', ir, lin, tadd, 2)
  plot_wavefun('sllleft', ir, lin, tadd, 3)
#"""
 

show()
