#!/usr/bin/env python

from numpy import *
from matplotlib.pyplot import *

t = loadtxt('tmat_atom_001_energ_001.dat')
t = t[:,0]+1j*t[:,1]
t = t.reshape(int(sqrt(len(t))), -1)
t = mat(t)

with open('tmat_atom_001_energ_001.dat') as f:
    tmpline = f.readlines()[0]
    etxt = tmpline.split('=')[-1].split()
    e = float(etxt[0])+1j*float(etxt[1])

figure(figsize=(6,8))
subplot(4,2,1); imshow(real(t)); colorbar(); title('Re'); ylabel('t')
subplot(4,2,2); imshow(imag(t)); colorbar(); title('Im')

dt = (t.transpose().conjugate() - t)/2.
subplot(4,2,3); imshow(real(dt)); colorbar(); ylabel('lhs=(t^dagger-t)/2')
subplot(4,2,4); imshow(imag(dt)); colorbar()

dt2 = 1j*sqrt(real(e))*(t.transpose().conjugate() * t)
subplot(4,2,5); imshow(real(dt2)); colorbar(); ylabel('rhs=i*sqrt(e)*t^dagger*t')
subplot(4,2,6); imshow(imag(dt2)); colorbar()

dt = dt-dt2
#dt[abs(dt2)<10**-2]=0
subplot(4,2,7); imshow(log10(abs(real(dt)))); colorbar(); ylabel('log10(abs(lhs-rhs))')
subplot(4,2,8); imshow(log10(abs(imag(dt)))); colorbar()

suptitle('Checks for optical theorem')
subplots_adjust(top=0.891, bottom=0.048, left=0.043, right=0.943, hspace=0.483, wspace=0.174)

t = array(t)
for il in range(len(t)):
  print il, log10(abs( imag(t[il,il])-sqrt(real(e))*sum(abs(t[:,il])**2) ))

show()
