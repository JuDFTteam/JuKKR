#!/usr/bin/env python

from numpy import loadtxt, sqrt, sum
from matplotlib.pyplot import imshow, colorbar, show, figure, subplot, axis, subplots_adjust

figure(figsize=(12,6))
subplots_adjust(left=0.04, bottom=0.25, right=0.99, top=0.75)

# first kind
subplot(1,2,1)
title('wronskian of first kind')
d = sum(loadtxt('test_wronskian')[:,2:], axis=1)
d = d.reshape(int(sqrt(len(d))),-1)
imshow(d, interpolation='nearest')
axis('equal')
colorbar()

# second kind
subplot(1,2,2)
title('wronskian of second kind')
d = sum(loadtxt('test_wronskian2')[:,2:], axis=1)
d = d.reshape(int(sqrt(len(d))),-1)
imshow(d, interpolation='nearest')
axis('equal')
colorbar()

show()
