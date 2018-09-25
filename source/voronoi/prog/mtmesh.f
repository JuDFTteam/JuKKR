      SUBROUTINE mtmesh(NRAD,npan,meshn,nm,xrn,drn,
     &                  nfu,thetas,lmifun,mtradius)
      Implicit None
c#@# KKRtags: VORONOI radial-grid
      include 'inc.geometry'
C
C
C Program  mtmesh.f adds one extra pannel inside the
C muffin-tin sphere to allow lattice relaxations.
C stores the mt-nized shapes in unit 15 as shapefun
C     .. Parameters ..
c     nrad : number of points added inside the MT radius
      Integer ibmaxd
      Parameter (IBMAXD=(LMAXD1+1)*(LMAXD1+1))
      Integer nrad
c      Parameter (nrad=20)
C     ..
C     .. Local Scalars ..
      REAL*8           dist,dn1,mtradius,pi,scale
      Integer ifun,ipan1,ir,lm,meshn,nfu,npan,number
      integer ip,i
C     ..
C     .. Arrays ..
      REAL*8           drn(IRID),thetas(IRID,ibmaxd),xrn(IRID)
      Integer nm(npand),lmIFUN(ibmaxd)
C     ..
C     .. Local Arrays ..      
      REAL*8           drn1(IRID),thetas1(IRID,ibmaxd),xrn1(IRID)
      Integer nm1(npand)
      integer meshn1,npan1
C     ..
C     .. Intrinsic Functions ..
      Intrinsic abs,datan,dsqrt,sqrt
C     ..
      pi = 4.d0*datan(1.d0)

      npan1 = npan + 1
      meshn1 = meshn + nrad
      nm1(1) = nrad
      
      do ip=2,npan1
         nm1(ip) = nm(ip-1)
      end do

      If (npan1.gt.npand) Then
        Write (6,FMT=*) ' npan , npand ',npan1,npand
        Stop
      End if
      
      If (meshn1.gt.irid) Then
        Write (6,FMT=*) ' meshn , irid ',meshn1,irid
        Stop
      End if

      dist = xrn(1) - mtradius

      If (dist.lt.1.0d-5) Then
         Write (6,FMT=*) 'Error from MTMESH '
         write (6,*) 
     & 'Your MT-radious is biger that the minimum shape radious ' 
         write(6,*) 'Your MT-Radious .....',mtradius
         write(6,*) 'Shape Radious .......',xrn(1)    
        Stop
      End if
           
      dn1 = dist/ (nrad-1)
      xrn1(nrad) = xrn(1)
      drn1(nrad) = dn1
      Do ir = 1,nrad-1
        xrn1(ir) = mtradius + dn1* (ir-1)
        drn1(ir) = dn1
      End do
      do i=1,meshn1-nrad
         xrn1(i+nrad) = xrn(i)
         drn1(i+nrad) = drn(i) 
      end do

      Do ir = 1,nrad
        thetas1(ir,1) = dsqrt(4.d0*pi)
        Do ifun = 2,ibmaxd
          thetas1(ir,ifun) = 0.0d0
        End do
      End do

      Do 10 ifun = 1,nfu
        do ir=1,meshn1-nrad
         thetas1(nrad+ir,ifun) = thetas(ir,ifun)
        end do
   10 Continue
c
c  Now map back and return. 
c
      npan = npan1
      meshn = meshn1

      do ip=1,npan
         nm(ip) = nm1(ip)
      end do  
      do ir=1,meshn
         xrn(ir) = xrn1(ir)
         drn(ir) = drn1(ir)
      end do

      Do ifun = 1,nfu
        do ir=1,meshn
         thetas(ir,ifun) = thetas1(ir,ifun)
        end do
      end do
c
c Now store on disk
c
c     write(15,FMT=9000) npan,meshn
c     Write (15,FMT=9000) (nm(ipan1),ipan1=1,npan)
c     Write (15,FMT=9010) (xrn(ir),drn(ir),ir=1,meshn)
c     Write (15,FMT=9000) nfu
c     Do ifun = 1,nfu
c        Write (15,FMT=9000) lmifun(ifun)
c        Write (15,FMT=9010) (thetas(ir,ifun),ir=1,meshn)
c     end do
      RETURN 
 9000 Format (16i5)
 9010 Format (4d20.12)
      End



