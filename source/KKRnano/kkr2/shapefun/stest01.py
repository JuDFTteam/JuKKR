"""Testing shapes"""

import voronoi08
import shapewrapper
import numpy as np

#Simple cubic testcase

a = np.array([1.,0.,0.])
b = np.array([0.,1.,0.])
c = np.array([0.,0.,1.])

# the if is to exclude the origin
rvec = np.array([n*a+m*b+k*c for n in range(-2,3) for m in range(-2,3) for k in range(-2,3) if abs(n)+abs(m)+abs(k) != 0],
                order='F').T

num = rvec.shape[1]

nfaced = num

print(nfaced)
weights = [1.] * nfaced

rmt,rout,volume,nface,a3,b3,c3,d3,nvert,xvert,yvert,zvert = voronoi08.voronoi08(num,rvec,30,1.,weights,1e-5,1e-5,False)

print("Muffin tin radius: {}".format(rmt))
print("Outer radius:      {}".format(rout))
print("Volume:            {}".format(volume))
print("Number of faces:   {}".format(nface))
print("Number of vertices {}".format(nvert[:nface]))

npoi = 125
nmin = 5
lmax = 4
dlt = 0.05
meshnd = 125
npand = 15
ibmaxd = 169

#rmt,rout,volume,nface,a3,b3,c3,d3,nvert,xvert,yvert,zvert = voronoi08(nvec,rvec,nvertmax,weight0,weight,tolvdist,tolarea,output,[nfaced])

#npan,nm,xrn,drn,meshn,thetas_s,lmifun_s,nfun = shapewrapper(npoi,aface,bface,cface,dface,nmin,nvertices,xvert,yvert,zvert,nface,lmax,dlt,ibmaxd,meshnd,npand,[nfaced,nvertd])
print nfaced
print nvert, len(nvert)
npan,nm,xrn,drn,meshn,thetas_s,lmifun_s,nfun = shapewrapper.shapewrapper(npoi,a3,b3,c3,d3,nmin,nvert,xvert,yvert,zvert,nface,lmax,dlt,ibmaxd,meshnd,npand)
