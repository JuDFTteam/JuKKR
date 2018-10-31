#!/usr/bin/env ipython
from pylab import *
from numpy import *
eps = 10.0**(-50.0)
print "This program plots the wave function R_Lambda1,Lambda2. It works only for lcut=3."
data=loadtxt("../../1Fe_0Cu_sratest2/fort.4000")
data2=loadtxt("../../1Fe_0Cu_sratest2/fort.4010")
rmesh=loadtxt("../../1Fe_0Cu_sratest2/fort.5000")
rll=data[:,::2]+1j*data[:,1::2]
rll2=data2[:,::2]+1j*data[:,1::2]
Lambda1   = int(raw_input("Lambda1="))
Lambda2   = int(raw_input("Lambda2="))
position  = (Lambda1-1)*32 + Lambda2-1
print "line of the file fort.4000 or fort.4010: "
print position
figure(1)
title('Big component')
plot(rmesh[:],real(rll[position,:]))
plot(rmesh[:],imag(rll[position,:]))
figure(2)
title('Small component')
plot(rmesh[:],real(rll2[position,:]))
plot(rmesh[:],imag(rll2[position,:]))

rmat=rll[:,400]
for i in range(size(rmat)):
    if(abs(rmat[i]) < eps):
        rmat[i] = 0.0
rmat=rmat.reshape([32,32])
matshow(log(abs(rmat))/log(10),vmin=-36,vmax=0)
title('Big component matrix for the 400th r value')
colorbar()

rmat2=rll2[:,400]
for i in range(size(rmat2)):
    if(abs(rmat2[i]) < eps):
        rmat2[i] = 0.0
rmat2=rmat2.reshape([32,32])
matshow(log(abs(rmat2))/log(10),vmin=-36,vmax=0)
title('Small component matrix for the 400th r value')
colorbar()


show()
quit()