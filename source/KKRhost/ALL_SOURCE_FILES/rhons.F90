SUBROUTINE rhons(den,df,drdi,gmat,ek,rho2ns,ipan,ircut,irmin, &    ! Added IRMIN 1.7.2014  &
        thetas,ifunm,lmsp,nsra,qns,pns,ar,cr,pz,fz,qz,  &
        sz,cleb,icleb,jend,iend,ekl,denlm,gflle_part)
!-----------------------------------------------------------------------

!     the charge density is developed in spherical harmonics :

!             rho(r) =   { rho(lm,r) * y(r,lm) }       (summed over lm)

!          rho(lm,r) =   { do rho(r) * y(r,lm)         (integrated over
!                                                           unit sphere)
!     in the case of spin-polarization :
!       the spin density is developed in spherical harmonics :

!            sden(r) =   { sden(lm,r) * y(r,lm) }      (summed over lm)

!         sden(lm,r) =   { do sden(r) * y(r,lm)        (integrated over
!                                                           unit sphere)
!     n(r,e) is developed in

!        n(r,e) = { y(r,l'm') * n(l'm',lm,r,e) * y(r,lm) }

!     therefore a faltung of n(l'm',lm,r,e) with the gaunt coeffients
!     has to be used to calculate the lm-contribution of the charge
!     density .


!     calculate the valence density of states , in the spin-polarized
!      case spin dependent .
!     recognize that the density of states is always complex also in
!      the case of "real-energy-integation" (ief>0) since in that case
!      the energy integration is done parallel to the real energy axis
!      but not on the real energy axis .
!     in the last energy-spin loop the l-contribution of the valence
!      charge is calculated .

!                               b.drittler   aug. 1988

!     modified for the use of shape functions

!     attention : irmin + 3 has to be less then imt
!                 if shape functions are used

!                               b.drittler   july 1989
!-----------------------------------------------------------------------
      use mod_DataTypes
      IMPLICIT NONE
!.. Parameters ..
      INCLUDE 'inc.p'
! *********************************************************************
! * For KREL = 1 (relativistic mode)                                  *
! *                                                                   *
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *********************************************************************

      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER LMPOTD,LMMAXD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      parameter (lmmaxd= (krel+1) * (lmaxd+1)**2)
      INTEGER LMAXD1
      PARAMETER (LMAXD1= LMAXD+1)
!..
!.. Scalar Arguments ..
      complex (kind=dp) DF,EK
      INTEGER IEND,IPAN,NSRA,IRMIN
!..
!.. Array Arguments ..
      complex (kind=dp) AR(LMMAXD,LMMAXD),CR(LMMAXD,LMMAXD),DEN(0:LMAXD1), &
                     EKL(0:LMAXD),FZ(IRMD,0:LMAXD),GMAT(LMMAXD,LMMAXD), &
                     PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),PZ(IRMD,0:LMAXD), &
                     QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),QZ(IRMD,0:LMAXD), &
                     SZ(IRMD,0:LMAXD) &
                    ,DENLM(LMMAXD)
#ifndef CPP_MPI
      complex (kind=dp) ENERG  ! lm-dos
#endif
      real (kind=dp) CLEB(*),DRDI(IRMD),RHO2NS(IRMD,LMPOTD), &
                       THETAS(IRID,NFUND)
      INTEGER ICLEB(NCLEB,4),IFUNM(*),IRCUT(0:IPAND), &
              JEND(LMPOTD,0:LMAXD,0:LMAXD),LMSP(*)
!..
!.. Local Scalars ..
      complex (kind=dp) DENNS,V1
      INTEGER IMT1,L,LM,M,IRMAX, LM1, LM2
!..
!.. Local Arrays ..
      complex (kind=dp) CDEN(IRMD,0:LMAXD),CDENNS(IRMD),EFAC(LMMAXD) &
                    ,CDENLM(IRMD,LMMAXD),CWR(IRMD,LMMAXD,LMMAXD)   &! lm-dos
                    ,GFLLE_PART(LMMAXD,LMMAXD)
!..
!.. External Functions ..
      LOGICAL OPT                          ! qdos
      EXTERNAL OPT                         ! qdos
!..
!.. External Subroutines ..
      EXTERNAL CSIMPK,RHOIN,RHOOUT
!..
!.. Intrinsic Functions ..
      INTRINSIC DBLE
!..
      real (kind=dp) PI
      PI = 4.D0*DATAN(1.D0)

!---> set up efac(lm) = sqrt(e))**l/(2l - 1)!!

efac(1) = 1.0D0
v1 = 1.0D0
DO  l = 1,lmaxd
  v1 = v1*ek/DBLE(2*l-1)
  DO  m = -l,l
    lm = l* (l+1) + m + 1
    efac(lm) = v1
  END DO
END DO

imt1 = ircut(1)
irmax = ircut(ipan)

CALL rhoout(cden,df,gmat,ek,pns,qns,rho2ns,thetas,ifunm,ipan,  &
    imt1,irmin,irmax,lmsp,cdenns,nsra,cleb,icleb,iend    &  ! Added IRMIN,IRMAX 1.7.2014  &
    ,cdenlm,cwr)  ! lm-dos

CALL rhoin(ar,cden,cr,df,gmat,ek,rho2ns,irmin,nsra,efac,pz,fz, &  ! Changed from irmind TO irmin 1.7.2014  &
    qz,sz,cleb,icleb,jend,iend,ekl  &
    ,cdenlm)  ! lm-dos  ! Attention, cwr does not go into rhoin, does lmlm-dos work properly?


!---> calculate complex density of states

DO  l = 0,lmaxd
  
!---> call integration subroutine
  
  CALL csimpk(cden(1,l),den(l),ipan,ircut,drdi)
END DO

DO  lm1 = 1,lmmaxd  ! lm-dos
  CALL csimpk(cdenlm(1,lm1),denlm(lm1),ipan,ircut,drdi)  ! lm-dos
  IF (opt('lmlm-dos').OR.opt('qdos    ').OR. &          ! lmlm-dos  &
        opt('LDA+U   ')) THEN                            ! LDAU
    DO  lm2 = 1,lmmaxd                               ! lmlm-dos
      CALL csimpk(cwr(1,lm1,lm2),gflle_part(lm1,lm2), & ! lmlm-dos  &
          ipan,ircut,drdi)                     ! lmlm-dos
    END DO
  endif                                                ! lmlm-dos
END DO

! Energy depends on EK and NSRA:                            ! lm-dos
!     IF (NSRA.EQ.1) EK = SQRT(E)                           ! lm-dos
!     IF (NSRA.EQ.2) EK = SQRT(E+E*E/ (CVLIGHT*CVLIGHT))    ! lm-dos
!     CVLIGHT=274.0720442D0                                 ! lm-dos
! Therefore the following is a good approximation           ! lm-dos
! for energies of a few Ryd:                                ! lm-dos
#ifndef CPP_MPI
IF (.NOT.opt('qdos    ')) THEN
  energ = ek**2                                         ! lm-dos
  WRITE(30,9000) real(energ, kind=dp),(-aimag(denlm(lm))/pi,lm=1,lmmaxd)
  9000 FORMAT(30E12.4)
endif  ! not qdos option
#endif


IF (ipan > 1) THEN
  CALL csimpk(cdenns,denns,ipan,ircut,drdi)
  den(lmaxd1) = denns
endif

RETURN
END SUBROUTINE rhons
