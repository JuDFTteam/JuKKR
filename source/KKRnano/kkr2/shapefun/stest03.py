import voronoi08
import numpy as np
import shapewrapper

#fcc cubic testcase

a = np.array([0.5,0.5,0.0])
b = np.array([0.5,0.0,0.5])
c = np.array([0.0,0.5,0.5])

# the if is to exclude the origin
rvec = np.array([n*a+m*b+k*c for n in range(-2,3) for m in range(-2,3) for k in range(-2,3) if abs(n)+abs(m)+abs(k) != 0],
                order='F').T

num = rvec.shape[1]

nfaced = num

print(nfaced)
weights = [1.] * nfaced

rvecnew = rvec.copy()

# set npand to sum of number of vertices of each face + number of faces
# set npoi to at least npand*nmin
# set meshnd to npoi

#npoi = 300    # increase from 125
nmin = 5
lmax = 4
dlt = 0.05
#npand = 100   # increase from 15
ibmaxd = 169
npoi = 125 # min. 125 points

#rmt,rout,volume,nface,a3,b3,c3,d3,nvert,xvert,yvert,zvert = voronoi08(nvec,rvec,nvertmax,weight0,weight,tolvdist,tolarea,output,[nfaced])

#npan,nm,xrn,drn,meshn,thetas_s,lmifun_s,nfun = shapewrapper(npoi,aface,bface,cface,dface,nmin,nvertices,xvert,yvert,zvert,nface,lmax,dlt,ibmaxd,meshnd,npand,[nfaced,nvertd])


def dotest(rvec):

    TOL = 1e-10

    rmt,rout,volume,nface,a3,b3,c3,d3,nvert,xvert,yvert,zvert = voronoi08.voronoi08(num,rvec,30,1.,weights,TOL,TOL,False)

    print("Muffin tin radius: {}".format(rmt))
    print("Outer radius:      {}".format(rout))
    print("Volume:            {}".format(volume))
    print("Number of faces:   {}".format(nface))
    print("Number of vertices {} Sum: {}".format(nvert[:nface], sum(nvert[:nface])))

    npand = sum(nvert[:nface]) + nface
    meshnd = npand*nmin

    npan,nm,xrn,drn,meshn,thetas_s,lmifun_s,nfun = shapewrapper.shapewrapper(npoi,a3,b3,c3,d3,nmin,nvert,xvert,yvert,zvert,nface,lmax,dlt,ibmaxd,meshnd,npand)
    
    print("Using {} points and maximum of {} panels, {} panels created.".format(meshn, npand, npan))
    print(npoi)
    return npan,nm,xrn,drn,meshn,thetas_s,lmifun_s,nfun


print '-'*79   
rvecnew = rvec + ((np.random.rand(*rvecnew.shape) - .5) / 2.)
npan,nm,xrn,drn,meshn,thetas_s,lmifun_s,nfun = dotest(rvecnew)
