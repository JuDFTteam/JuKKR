import voronoi08
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

r = voronoi08.voronoi08(num,rvec,num,1.,weights,1e-5,1e-5,False)

print("Muffin tin radius: {}".format(r[0]))
print("Outer radius:      {}".format(r[1]))
print("Volume:            {}".format(r[2]))
print("Number of faces:   {}".format(r[3]))
print("Number of vertices {}".format(r[8][:r[3]]))

