subroutine rhons(den,df,drdi,gmat,ek,rho2ns,ipan,ircut,thetas, &
ifunm,lmsp,nsra,qns,pns,ar,cr,pz,fz,qz,sz,cleb, &
icleb,jend,iend,ekl, &
lmax, irmd, irnsd, irid, ipand, nfund, ncleb)
  !-----------------------------------------------------------------------
  !
  !     the charge density is developed in spherical harmonics :
  !
  !             rho(r) =   { rho(lm,r) * y(r,lm) }       (summed over lm)
  !
  !          rho(lm,r) =   { do rho(r) * y(r,lm)         (integrated over
  !                                                           unit sphere)
  !     in the case of spin-polarization :
  !       the spin density is developed in spherical harmonics :
  !
  !            sden(r) =   { sden(lm,r) * y(r,lm) }      (summed over lm)
  !
  !         sden(lm,r) =   { do sden(r) * y(r,lm)        (integrated over
  !                                                           unit sphere)
  !     n(r,e) is developed in
  !
  !        n(r,e) = { y(r,l'm') * n(l'm',lm,r,e) * y(r,lm) }
  !
  !     therefore a faltung of n(l'm',lm,r,e) with the gaunt coeffients
  !     has to be used to calculate the lm-contribution of the charge
  !     density .
  !
  !
  !     calculate the valence density of states , in the spin-polarized
  !      case spin dependent .
  !     recognize that the density of states is always complex also in
  !      the case of "real-energy-integation" (ief>0) since in that case
  !      the energy integration is done parallel to the real energy axis
  !      but not on the real energy axis .
  !     in the last energy-spin loop the l-contribution of the valence
  !      charge is calculated .
  !
  !                               b.drittler   aug. 1988
  !
  !     modified for the use of shape functions
  !
  !     attention : irmin + 3 has to be less then imt
  !                 if shape functions are used
  !
  !                               b.drittler   july 1989
  !-----------------------------------------------------------------------
  implicit none

  integer lmax
  integer irmd
  integer ncleb
  integer irnsd
  integer irid
  integer ipand
  integer nfund

  !     INTEGER IRMIND
  !     PARAMETER (IRMIND=IRMD-IRNSD)
  !     INTEGER LMPOTD,LMMAXD
  !     PARAMETER (LMPOTD= (LPOTD+1)**2) ! LMPOTD= (2*LMAX+1)**2
  !     PARAMETER (LMMAXD= (LMAXD+1)**2)
  !     INTEGER LMAXD1
  !     PARAMETER (LMAXD1= LMAXD+1)
  !     ..
  !     .. Scalar Arguments ..
  double complex df,ek
  integer iend,ipan,nsra
  !     ..
  !     .. Array Arguments ..
  !     DOUBLE COMPLEX AR(LMMAXD,LMMAXD),CR(LMMAXD,LMMAXD),DEN(0:LMAXD1),
  !    +               EKL(0:LMAXD),FZ(IRMD,0:LMAXD),GMAT(LMMAXD,LMMAXD),
  !    +               PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),PZ(IRMD,0:LMAXD),
  !    +               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),QZ(IRMD,0:LMAXD),
  !    +               SZ(IRMD,0:LMAXD)
  !     DOUBLE PRECISION CLEB(*),DRDI(IRMD),RHO2NS(IRMD,LMPOTD),
  !    +                 THETAS(IRID,NFUND)
  !     INTEGER ICLEB(NCLEB,3),IFUNM(*),IRCUT(0:IPAND),
  !    +        JEND(LMPOTD,0:LMAXD,0:LMAXD),LMSP(*)

  double complex ar((lmax+1)**2,(lmax+1)**2)
  double complex cr((lmax+1)**2,(lmax+1)**2)
  double complex den(0:lmax+1)
  double complex ekl(0:lmax)
  double complex fz(irmd,0:lmax)
  double complex gmat((lmax+1)**2,(lmax+1)**2)
  double complex pns((lmax+1)**2,(lmax+1)**2,irmd-irnsd:irmd,2)
  double complex pz(irmd,0:lmax)
  double complex qns((lmax+1)**2,(lmax+1)**2,irmd-irnsd:irmd,2)
  double complex qz(irmd,0:lmax)
  double complex sz(irmd,0:lmax)

  !     DOUBLE PRECISION CLEB(*),DRDI(IRMD),RHO2NS(IRMD,LMPOTD),
  !    +                 THETAS(IRID,NFUND)
  !     INTEGER ICLEB(NCLEB,3),IFUNM(*),IRCUT(0:IPAND),
  !    +        JEND(LMPOTD,0:LMAXD,0:LMAXD),LMSP(*)

  double precision cleb(*)
  double precision drdi(irmd)
  double precision rho2ns(irmd,(2*lmax+1)**2)
  double precision thetas(irid,nfund)

  integer icleb(ncleb,3)
  integer ifunm(*)
  integer ircut(0:ipand)
  integer jend((2*lmax+1)**2,0:lmax,0:lmax)
  integer lmsp(*)
  !     ..
  !     .. Local Scalars ..
  double complex denns,v1
  integer imt1,l,lm,m
  !     ..
  !     .. Local Arrays ..
  double complex cden(irmd,0:lmax)
  double complex cdenns(irmd)
  double complex efac((lmax+1)**2)
  !     ..
  !     .. External Subroutines ..
  external csimpk,rhoin,rhoout
  !     ..
  !     .. Intrinsic Functions ..
  intrinsic dble

  integer irmind
  irmind=irmd-irnsd
  !
  !---> set up efac(lm) = sqrt(e))**l/(2l - 1)!!
  !
  efac(1) = 1.0d0
  v1 = 1.0d0
  do 20 l = 1,lmax
    v1 = v1*ek/dble(2*l-1)
    do 10 m = -l,l
      lm = l* (l+1) + m + 1
      efac(lm) = v1
10  continue
20 continue
   !
   imt1 = ircut(1)
   !

   call rhoout(cden,df,gmat,ek,pns,qns,rho2ns,thetas,ifunm,ipan, &
   imt1,lmsp,cdenns,nsra,cleb,icleb,iend, &
   lmax, irmd, irnsd, irid, nfund, ncleb)

   !
   call rhoin(ar,cden,cr,df,gmat,ek,rho2ns,irmind,nsra,efac,pz,fz, &
   qz,sz,cleb,icleb,jend,iend,ekl, &
   lmax, irmd, ncleb)

   !
   !---> calculate complex density of states
   !
   do 30 l = 0,lmax
     !
     !---> call integration subroutine
     !
     call csimpk(cden(1,l),den(l),ipan,ircut,drdi)
30 continue

   if (ipan.gt.1) then
     call csimpk(cdenns,denns,ipan,ircut,drdi)
     den(lmax+1) = denns
   end if

   return
 end
