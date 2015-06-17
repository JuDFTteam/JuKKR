!>    Convolute the potential with the shapefunctions.
!>    Operates on one spin channel and one site.
!>    @param GSH shape-Gaunt-coefficients
!>    @param Z    nuclear charge
!>    @param R    radial mesh
!>    @param VONS non-spherical part of potential (only for one spin
!>                channel)
! **********************************************************************
subroutine convol_new(imt1,irc1,imaxsh,ilm,ifunm,lmpot,gsh, &
thetas,z,r,vons,lmsp, &
irid, nfund, irmd, ngshd)
! **********************************************************************
  implicit none

  !      INTEGER LMPOTD
  !      PARAMETER (LMPOTD= (LPOTD+1)**2)

  !     Number of radial mesh points for shape functions
  integer irid
  !     maximal number of shape functions
  integer nfund
  integer irmd
  integer ngshd
  !     ..
  !     .. Scalar Arguments ..
  double precision rfpi,z
  integer imaxsh,imt1,irc1,lmpot
  !     ..
  !     .. Array Arguments ..
  double precision gsh(*),r(*),thetas(irid,nfund),vons(irmd,*)
  integer ifunm(*),ilm(ngshd,3),lmsp(*)
  !     ..
  !     .. Local Scalars ..
  double precision zzor
  integer i,ifun,ir,irh,lm,lm1,lm2,lm3
  !     ..
  !     .. Local Arrays ..
  double precision vstore(irid,lmpot)
  !     ..

  rfpi = sqrt(16.0d0 * atan(1.0d0))

  do 20 lm = 1,lmpot
    do 10 ir = 1,irc1 - imt1
      vstore(ir,lm) = 0.0d0
10  continue
20 continue

   do 30 ir = imt1 + 1,irc1
     zzor = 2.0d0*z/r(ir)*rfpi
     vons(ir,1) = vons(ir,1) - zzor
30 continue

   do 50 i = 1,imaxsh
     lm1 = ilm(i,1)
     lm2 = ilm(i,2)
     lm3 = ilm(i,3)
     if (lmsp(lm3).gt.0) then
       ifun = ifunm(lm3)
       do 40 ir = imt1 + 1,irc1
         irh = ir - imt1
         vstore(irh,lm1) = vstore(irh,lm1) + &
         gsh(i)*vons(ir,lm2)*thetas(irh,ifun)
40     continue
     end if
50 continue

   do 60 ir = imt1 + 1,irc1
     irh = ir - imt1
     zzor = 2.0d0*z/r(ir)*rfpi
     vons(ir,1) = vstore(irh,1) + zzor
60 continue

   do 80 lm = 2,lmpot
     do 70 ir = imt1 + 1,irc1
       irh = ir - imt1
       vons(ir,lm) = vstore(irh,lm)
70   continue
80 continue

 end
