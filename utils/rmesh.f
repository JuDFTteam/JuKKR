      program radmesh
c Gives the exponential mesh used in the KKR program.
c Formula:
c r(i) = b * (exp(a*(i-1))-1)
c r(0) = 0, r(IMT) = RMT
c
c Note that the first non-zero point is at
c
c rmesh(2) = b*(exp(a * 2.)-1.)
c with 
c b=rmt/(exp(a*(imt-1))-1)
c i.e.
c rmesh(2) = rmt * (exp(a * 2.)-1.) / (exp(a*(imt-1))-1)
c
      implicit none
c Parameters/dimensions:
      integer ndim
      parameter (ndim=10000)
c Input
      integer imt  ! Mesh endpoint
c Input/Output (depending on mode)
      real*8 rmt,aa,bb,r2
c Output
      real*8 rmesh(ndim)
c Local
      integer ipoint

 100  write(*,*) 'Input Imt,Rmt,a,b'
      write(*,*) '(One of Rmt,a,b is derivable from the others,'
      write(*,*) ' input for this should be zero;'
      write(*,*) ' the other two should be positive.)'
      
      read(*,*) imt,rmt,aa,bb

      if (rmt.ne.0.d0.and.aa.ne.0.d0.and.bb.ne.0.d0) then
         write(*,*) 'Wrong input: Rmt,a,b>0',imt,rmt,aa,bb
         goto 100
      endif

      if (imt.le.1.or.imt.gt.ndim) then 
         write(*,*) 'Wrong input: Imt ',imt
         goto 100
      endif

      if (aa.eq.0.d0) then
         aa = dlog(rmt / bb + 1.d0) / dfloat(imt - 1)
      else if (bb.eq.0.d0) then
         bb = rmt / (dexp(aa*dfloat(imt-1))-1)
      else
         rmt = bb * (dexp(aa*dfloat(imt-1))-1)
      endif

      do ipoint = 1,imt
         rmesh(ipoint) = bb * (dexp(aa*dfloat(ipoint-1))-1)
      enddo

      write(1,8000) '# Imt,Rmt,a,b: ',imt,rmt,aa,bb
      do ipoint = 1,imt
         write(1,9000) ipoint,rmesh(ipoint)
      enddo

 8000 format(a15,i5,3e16.8)
 9000 format(i5,e16.8)

      end
