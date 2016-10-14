from numpy import *
import sys,mynormalize
from matplotlib.pyplot import figure,subplot,pcolormesh,xlim,ylim,colorbar,title,show,gci,clim,axis

path = 'spinpol_Bi2Te3/SYSTEM2/'

out = [0,0,0,0,0,0]
for i in range(len(out)):
            out[i] = load(path+'saved_data_%i.npy'%i)
x,y,dat,sx,sy,sz = out
dat,sx,sy,sz = transpose(dat),transpose(sx),transpose(sy),transpose(sz)


c_scale = 1.0

d = zeros((300,300))
dsx= zeros((300,300))
dsy= zeros((300,300))
dsz= zeros((300,300))
d[100:200,100:200] = dat
dsx[100:200,100:200] = sx
dsy[100:200,100:200] = sy
dsz[100:200,100:200] = sz
dat = d
sz  = dsz
sy  = dsy
sx  = dsx
x = linspace(3*x.min(),3*x.max(),len(A))
y = linspace(3*y.min(),3*y.max(),len(A))



#'''
figure()
subplot(2,2,1)
title('log(I)')
pcolormesh(x,y,log(dat))
xlim(x.min(),x.max())
ylim(y.min(),y.max())
colorbar()
subplot(2,2,2)
title('sx')
pcolormesh(x,y,sx,cmap='seismic')
xlim(x.min(),x.max())
ylim(y.min(),y.max())
colorbar()
subplot(2,2,3)
title('sy')
pcolormesh(x,y,sy,cmap='seismic')
xlim(x.min(),x.max())
ylim(y.min(),y.max())
colorbar()
subplot(2,2,4)
title('sz')
pcolormesh(x,y,sz,cmap='seismic')
xlim(x.min(),x.max())
ylim(y.min(),y.max())
colorbar()
#'''

n = sqrt(sx**2+sy**2+sz**2)
n[n==0] = 1
sx,sy,sz = sx/n,sy/n,sz/n

def self_convol(dat):
        return real(fft.fftshift(fft.ifft2(conjugate(fft.fft2(dat))*fft.fft2(dat))))

def get_phi(r):       
      r0 = [1,0]
      sg = sign(cross([r0[0],r0[1],0],[r[0],r[1],0])[-1])
      phi = sg*nan_to_num(arccos((r0[0]*r[0]+r0[1]*r[1])/(sqrt(r[0]**2+r[1]**2)*sqrt(r0[0]**2+r0[1]**2))))
      return phi


def xy_proj(dat):
        l = shape(dat)
        outx = zeros_like(dat)
        outy = zeros_like(dat)
        for i in xrange(l[0]):
                for j in xrange(l[1]):
                        r_ij = [i-l[0]/2,j-l[1]/2]
                        outx[i,j] = cos(get_phi(r_ij))
                        outy[i,j] = sin(get_phi(r_ij))
        outx[:l[0]/2,l[1]/2]-=2
        return outx,outy

figure()
title('jdos(q)')
jdos=transpose(self_convol(dat))
jdos = jdos/max(jdos.reshape(-1))
jdos[jdos==1]=0
pcolormesh(jdos)
colorbar()

figure()
title('spp(q)')
xmat,ymat = xy_proj(dat)
A  = transpose(self_convol(ymat*dat)+self_convol(xmat*dat)+self_convol(dat))
A = A/max(A.reshape(-1))
A[A==1] = 0
img = pcolormesh(x,y,A)
cbar = colorbar()
axis('equal')
xlim(-0.35,0.35);ylim(-0.35,0.35)
show()

"""
figure()
title('jdos(q) - spp(q)')
img = pcolormesh(A-jdos)
cbar = colorbar()
show()
"""


"""
subplot(2,2,1)
title('I(k)')
pcolormesh(dat)
colorbar()
"""
"""
figure()
subplot(2,2,1)
title('jdos(q)=selfconvol(I)')
jdos=self_convol(dat)
pcolormesh(jdos)
colorbar()
clims = gci().get_clim()
print clims[0],clims[1]*c_scale
clims = real(clims)
print clims[0],clims[1]*c_scale
print clims[0]<clims[1]*c_scale
clim(clims[0],clims[1]*c_scale)
subplot(2,2,2)
title('log(jdos)')
jdos[jdos==0]=10**-10
jdos=log(jdos)
pcolormesh(jdos)
colorbar()

xmat,ymat = xy_proj(dat)
'''
subplot(2,2,3)
title('x projection')
pcolormesh(xmat)
colorbar()
subplot(2,2,4)
title('y projection')
pcolormesh(ymat)
colorbar()


figure()
subplot(2,2,1)
title('I(k)')
pcolormesh(transpose(dat))
colorbar()
xmat,ymat = xy_proj(dat)
subplot(2,2,2)
title('xmat*I(k)')
A  = array(xmat*dat)
pcolormesh(A)
colorbar()
subplot(2,2,3)
title('ymat*I(k)')
A  = array(ymat*dat)
pcolormesh(A)
colorbar()
'''
#subplot(2,2,4)
subplot(2,2,3)
title('spp = selfconvol(ymat*I)+selfconvol(xmat*I)+selfconvol(I)')
A  = self_convol(ymat*dat)+self_convol(xmat*dat)+self_convol(dat)
img = pcolormesh(A)
cbar = colorbar()
clims = gci().get_clim()
clim(clims[0],clims[1]*c_scale)
subplot(2,2,4)
title('log(spp)')
A[A==0] = 10**-8
pcolormesh(log(A))
colorbar()


figure()
subplot(2,2,1)
title('spp')
pcolormesh(A)
colorbar()
subplot(2,2,2)
title('spp''_xy = selfconvol(sx*I)+selfconvol(sy*I)')
A  = self_convol(sx*dat)+self_convol(sy*dat)#+self_convol(sz*dat)#+self_convol(dat)
pcolormesh(A)
colorbar()
clims = gci().get_clim()
clim(clims[0],clims[1]*c_scale)
subplot(2,2,3)
title('spp''_xyz = selfconvol(sx*I)+selfconvol(sy*I)+selfconvol(sz*I)')
A  = self_convol(sx*dat)+self_convol(sy*dat)+self_convol(sz*dat)#+self_convol(dat)
img = pcolormesh(A)
cbar = colorbar()
clims = gci().get_clim()
clim(clims[0],clims[1]*c_scale)
#cbar.set_norm(mynormalize.MyNormalize(vmin=A.min(),vmax=A.max(),stretch='linear'))
#cbar = mynormalize.DraggableColorbar(cbar,img)
#cbar.connect()
subplot(2,2,4)
title('spp'' = selfconvol(sx*I)+selfconvol(sy*I)+selfconvol(sz*I)+selfconvol(I)')
A  = self_convol(sx*dat)+self_convol(sy*dat)+self_convol(sz*dat)+self_convol(dat)
pcolormesh(A)
colorbar()
clims = gci().get_clim()
clim(clims[0],clims[1]*c_scale)



figure()
subplot(2,2,1)
title('spp'' xpart')
A = self_convol(sx*dat)
pcolormesh(A)
colorbar()
clims = gci().get_clim()
clim(clims[0],clims[1]*c_scale)
subplot(2,2,2)
title('spp'' ypart')
A  = self_convol(sy*dat)#+self_convol(sz*dat)#+self_convol(dat)
pcolormesh(A)
colorbar()
clims = gci().get_clim()
clim(clims[0],clims[1]*c_scale)
subplot(2,2,3)
title('spp'' zpart')
A  = self_convol(sz*dat)#+self_convol(dat)
pcolormesh(A)
colorbar()
clims = gci().get_clim()
clim(clims[0],clims[1]*c_scale)
subplot(2,2,4)
title('spp'' jdos-part')
A  = self_convol(dat)
pcolormesh(A)
colorbar()
clims = gci().get_clim()
clim(clims[0],clims[1]*c_scale)



figure()
subplot(2,2,1)
title('spp xpart')
A = self_convol(xmat*dat)
pcolormesh(A)
colorbar()
clims = gci().get_clim()
clim(clims[0],clims[1]*c_scale)
subplot(2,2,2)
title('spp ypart')
A  = self_convol(ymat*dat)#+self_convol(sz*dat)#+self_convol(dat)
pcolormesh(A)
colorbar()
clims = gci().get_clim()
clim(clims[0],clims[1]*c_scale)
subplot(2,2,3)
title('spp jdos-part')
A  = self_convol(dat)
pcolormesh(A)
colorbar()
clims = gci().get_clim()
clim(clims[0],clims[1]*c_scale)
subplot(2,2,4)
title('I(k)')
A  = dat
pcolormesh(A)
colorbar()


'''
figure()
subplot(2,2,1)
title('log(ifft(spp''))')
pcolormesh(log(abs(fft.fftshift(fft.ifft2(A)))))
colorbar()
subplot(2,2,2)
title('real(ifft(spp''))')
pcolormesh(real(fft.fftshift(fft.ifft2(A))))
colorbar()
subplot(2,2,3)
title('imag(ifft(spp''))')
pcolormesh(imag(fft.fftshift(fft.ifft2(A))))
colorbar()
subplot(2,2,4)
title('ifft(spp'')')
pcolormesh(fft.fftshift(fft.ifft2(A)))
colorbar()
'''




show()
"""
