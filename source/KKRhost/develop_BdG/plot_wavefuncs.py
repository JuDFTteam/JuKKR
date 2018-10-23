#!/usr/bin/env ipython

from numpy import *
from matplotlib.pyplot import *

# source terms
#"""
figure()
dim  = [int(i.split()[0]) for i in array(open('rll_source_jlk_atom_001_energ_001.dat').readline().split('='))[[1,2]]]

subplot(2,2,1)
jlk = loadtxt('rll_source_jlk_atom_001_energ_001.dat')
jlk = (jlk[:,0]+1j*jlk[:,1]).reshape(dim[1], dim[0])
plot(abs(jlk))
title('abs(jlk)')

subplot(2,2,2)
hlk = loadtxt('rll_source_hlk_atom_001_energ_001.dat')
hlk = (hlk[:,0]+1j*hlk[:,1]).reshape(dim[1], dim[0])
plot(abs(hlk))
title('abs(hlk)')

subplot(2,2,3)
jlk2= loadtxt('rll_source_jlk2_atom_001_energ_001.dat')
jlk2= (jlk2[:,0]+1j*jlk2[:,1]).reshape(dim[1], dim[0])
plot(abs(jlk2))
title('abs(jlk2)')

subplot(2,2,4)
hlk2= loadtxt('rll_source_hlk2_atom_001_energ_001.dat')
hlk2= (hlk2[:,0]+1j*hlk2[:,1]).reshape(dim[1], dim[0])
plot(abs(hlk))
title('abs(hlk2)')
#"""


# rll/sll
dim  = [int(i.split()[0]) for i in array(open('rll_atom_001_energ_001.dat').readline().split('='))[[1,2,3]]]

rll = loadtxt('rll_atom_001_energ_001.dat')
sll = loadtxt('sll_atom_001_energ_001.dat')
rll = (rll[:,0]+1j*rll[:,1]).reshape(dim[2], dim[1], dim[0])
sll = (sll[:,0]+1j*sll[:,1]).reshape(dim[2], dim[1], dim[0])

# take trace in lm (big component only)
out = []
for j in range(len(rll)):
  tmpsum = 0+1j*0
  for i in range(len(rll[0])):
    tmpsum+= rll[j,i,i]
  out.append(tmpsum)
rll_tr = array(out)
out = []
for j in range(len(rll)):
  tmpsum = 0+1j*0
  for i in range(len(rll[0])):
    tmpsum+= sll[j,i,i]
  out.append(tmpsum)
sll_tr = array(out)

figure()
subplot(2,2,1)
plot(abs(rll_tr))
plot(abs(sll_tr), '--')
title('abs(Tr[rll])')
subplot(2,2,3)
plot(real(rll_tr))
plot(real(sll_tr), '--')
title('real(Tr[rll])')
subplot(2,2,4)
plot(imag(rll_tr))
plot(imag(sll_tr), '--')
title('imag(Tr[rll])')
suptitle('big component')

# take trace in lm (big component only)
out = []
for j in range(len(rll)):
  tmpsum = 0+1j*0
  for i in range(len(rll[0])):
    tmpsum+= rll[j,i,i+dim[1]]
  out.append(tmpsum)
rll_tr = array(out)
out = []
for j in range(len(rll)):
  tmpsum = 0+1j*0
  for i in range(len(rll[0])):
    tmpsum+= sll[j,i,i+dim[1]]
  out.append(tmpsum)
sll_tr = array(out)

figure()
subplot(2,2,1)
plot(abs(rll_tr))
plot(abs(sll_tr), '--')
title('abs(Tr[rll])')
subplot(2,2,3)
plot(real(rll_tr))
plot(real(sll_tr), '--')
title('real(Tr[rll])')
subplot(2,2,4)
plot(imag(rll_tr))
plot(imag(sll_tr), '--')
title('imag(Tr[rll])')
suptitle('small component')

# some l-diagonal blocks
out = []
for j in range(len(rll)):
  for i in range(len(rll[0])):
    out.append(rll[j,i,i])
rll_ii = array(out).reshape(len(rll), len(rll[0]))
out = []
for j in range(len(rll)):
  for i in range(len(rll[0])):
    out.append(sll[j,i,i])
sll_ii = array(out).reshape(len(rll), len(rll[0]))

figure()
subplot(2,2,1)
for i in range(9):
   plot(abs(rll_ii[:,i]), label=i)
   plot(abs(sll_ii[:,i]), '--', label=i)
legend(fontsize='x-small')
title('abs(rll[i,i])')

subplot(2,2,3)
for i in range(9):
   plot(real(rll_ii[:,i]), label=i)
   plot(real(sll_ii[:,i]), '--', label=i)
legend(fontsize='x-small')
title('real(rll[i,i])')

subplot(2,2,4)
for i in range(9):
   plot(imag(rll_ii[:,i]), label=i)
   plot(imag(sll_ii[:,i]), '--', label=i)
legend(fontsize='x-small')
title('imag(rll[i,i])')

# sum over radial dimension:
rll_int = (sum(abs(rll), axis=0))
sll_int = (sum(abs(sll), axis=0))
figure()
subplot(2,1,1)
imshow(rll_int)
colorbar()
subplot(2,1,2)
imshow(sll_int)
colorbar()

rll_int = log(sum(abs(rll), axis=0))
sll_int = log(sum(abs(sll), axis=0))
figure()
subplot(2,1,1)
imshow(rll_int)
colorbar()
subplot(2,1,2)
imshow(sll_int)
colorbar()


show()
