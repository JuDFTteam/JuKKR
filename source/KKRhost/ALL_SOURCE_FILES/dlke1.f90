! 01.06.99 *************************************************************
SUBROUTINE dlke1(gllke,alat,nacls,naclsmax,rr,ezoa, atom,bzkp,ic,ginp,rcls)
! **********************************************************************

!     Fourier transformation of the cluster Greens function GINP

! ----------------------------------------------------------------------
use mod_types, only: t_inc
implicit none
!.. Parameters ..
include 'inc.p'
!
! *********************************************************************
! * For KREL = 1 (relativistic mode)                                  *
! *                                                                   *
! *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     *
! *          function, set up in the spin-independent non-relativstic *
! *          (l,m_l)-representation                                   *
! *                                                                   *
! *********************************************************************
!
INTEGER LMGF0D
PARAMETER (LMGF0D= (LMAXD+1)**2)
INTEGER ALMGF0
PARAMETER (ALMGF0=LMGF0D*NAEZD)
DOUBLE COMPLEX CI
PARAMETER (CI= (0.0D0,1.0D0))
!..
!.. Scalar Arguments ..
      DOUBLE PRECISION ALAT
      INTEGER IC,NACLSMAX
!..
!.. Array Arguments ..
INTEGER ATOM(*), &
        EZOA(*), &
        NACLS(*)
DOUBLE COMPLEX GLLKE(ALMGF0,*), &
               GINP(LMGF0D*NACLSMAX,*)
DOUBLE PRECISION BZKP(*), &
                 RR(3,0:NRD), &
                 RCLS(3,*)
!..
!.. Local Scalars ..
      DOUBLE PRECISION CONVPU,TPI
      INTEGER AM,I,II,IM,LM2,M
      DOUBLE COMPLEX  EIKR,TT
      LOGICAL OPT,TEST
!..
!.. Local Arrays ..
      DOUBLE COMPLEX ARG(3)
!..
!.. External Subroutines ..
      EXTERNAL CINIT,TEST,OPT,ZAXPY
!..
!.. Intrinsic Functions ..
      INTRINSIC ATAN,EXP
!..
!.. Save statement ..
      SAVE
!..
ii = 3
IF (opt('COMPLEX ')) ii = 6
IF (test('BZKP    ').AND.(t_inc%i_write>0))  &
    WRITE(1337,FMT='(6f12.6)') (bzkp(i),i=1,ii)

tpi = 8.0D0*ATAN(1.0D0)                 ! = 2*PI
convpu = alat/tpi

CALL cinit(lmgf0d*naezd*lmgf0d,gllke)

DO  m = 1,nacls(ic)
  
  
! --->  for option 'WIRE': avoid artifical couplings in the structural
!       Greens Function in in-plane-direction (perp. to c-axis)
  
  IF (atom(m) < 0) CYCLE
  
  IF (opt('ONEBULK ')) THEN     ! added 1.02.2000
!                                       corrected on 25.02.2000
!     if the phase factor exp(ik(r-r')) is included      ~
!     in the G...so if we resolve the dyson part for the G
!     and not for the G (see Peter Lang Ph.D thesis)
    
!     Here we do   --                           nn'
!                  \                            ii'          ii'
!                  /  exp( ik(X  -X + R  -R  ))G   (E)  =   G   (k,E)
!                  --          n'  n   i'  i    LL'          LL'
!                  n'
    
! In this case rcls is always (by constraction symmetric around each
! atom this means that a minus sign will not affect the result of the
! summation
    
    arg(1) = -ci*tpi*rcls(1,m)
    arg(2) = -ci*tpi*rcls(2,m)
    arg(3) = -ci*tpi*rcls(3,m)
  ELSE
    
!     Here we do   --                  nn'
!                  \                   ii'          ii'
!                  /  exp(+ik(x  -x ))G   (E)  =   G   (k,E)
!                  --          n'  n   LL'          LL'
!                  n'
!  Be carefull a minus sign must be included here. RR is not
!  symmetric around each atom. The minus comes from the fact that
!  the repulsive potential GF is calculated for 0n and not n0!
!  and that is why we nead a minus sign extra!
    
    
    arg(1) = -ci*tpi*rr(1,ezoa(m))
    arg(2) = -ci*tpi*rr(2,ezoa(m))
    arg(3) = -ci*tpi*rr(3,ezoa(m))
  END IF
  
  tt = bzkp(1)*arg(1)+bzkp(2)*arg(2)+bzkp(3)*arg(3)
  
  IF (opt('COMPLEX ')) THEN
    tt = tt + ci*(bzkp(4)*arg(1)+bzkp(5)*arg(2)+bzkp(6)*arg(3))
  END IF
  
!        write(6,*) BZKP(1),ARG(1),BZKP(2),ARG(2),BZKP(3),ARG(3)
!        write(6,*) BZKP(4),BZKP(5),BZKP(6)
!        write(6,*) 'm,atom(m),tt',m,atom(m),tt
  
  eikr = EXP(tt) * convpu    ! convert to p.u.
  
!        write(6,*) 'eikr',eikr
  
  im = 1 + (m-1)      *lmgf0d
  am = 1 + (atom(m)-1)*lmgf0d
  DO  lm2 = 1,lmgf0d
    CALL zaxpy(lmgf0d,eikr,ginp(im,lm2),1,gllke(am,lm2),1)
  END DO
  
END DO

RETURN
9000 FORMAT(3F12.4)
9010 FORMAT(2F18.10)
END SUBROUTINE dlke1
