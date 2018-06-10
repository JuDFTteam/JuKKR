SUBROUTINE corehff(kap1,kap2,mj,s,nsol,bhf,gck,fck,rc,drdic,rnuc,  &
        nzero,nrc)
!   ********************************************************************
!   *                                                                  *
!   *   CALCULATE THE RELATIVISTIC HYPERFINEFIELDS FOR THE             *
!   *                  CURRENT  CORE STATE S                           *
!   *                                                                  *
!   *   THE WAVE FUNCTION  {G(K,S),F(K,S)}  IS NORMALIZED TO 1         *
!   *                                                                  *
!   ********************************************************************

IMPLICIT NONE


! PARAMETER definitions
REAL*8 E0,A0,CAUTOG

!CONVERSION FACTOR FOR HYPERFINE FIELDS FROM A.U. TO GAUSS
!                                 ELECTRON CHARGE     IN ESU
!                                 BOHR-RADIUS         IN CM


PARAMETER (E0=1.6021892D-19*2.997930D+09,A0=0.52917706D-08, &
           CAUTOG=E0/(A0*A0))

! Dummy arguments
INTEGER KAP1,KAP2,NRC,NSOL,NZERO,S
REAL*8 MJ,RNUC
REAL*8 BHF(2,2),DRDIC(NRC),FCK(2,2,NRC),GCK(2,2,NRC),RC(NRC)

! Local variables
REAL*8 AME(2,2),XX(5),YI(NRC),YY(5),ZI(NRC)
DOUBLE PRECISION DBLE,DSQRT
INTEGER I,K1,K2,N
REAL*8 YLAG

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ANGULAR HYPERFINE MATRIX ELEMENTS   SEE E.G.  E.M.ROSE
!        THE FACTOR  I  HAS BEEN OMITTED

ame(1,1) = 4.0D0*kap1*mj/(4.0D0*kap1*kap1-1.0D0)
IF ( nsol == 2 ) THEN
  ame(2,2) = 4.0D0*kap2*mj/(4.0D0*kap2*kap2-1.0D0)
  ame(1,2) = DSQRT(0.25D0-(mj/DBLE(kap1-kap2))**2)
  ame(2,1) = ame(1,2)
END IF
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

DO k1 = 1,nsol
  DO k2 = 1,nsol
    DO n = 1,nzero
      yi(n) = drdic(n) *(gck(k1,s,n)*fck(k2,s,n)+fck(k1,s,n)*gck(k2,s,n)  &
          )
    END DO
    CALL rint4pts(yi,nzero,zi)
    IF ( rnuc /= 0.0D0 ) THEN
      DO i = 1,5
        xx(i) = rc(nzero-5+i)
        yy(i) = zi(nzero-5+i)
      END DO
      zi(nzero) = ylag(rnuc,xx,yy,0,4,5)
    END IF
    xx(1) = 1.0D0
!                      !RC( 1)
    xx(2) = 6.0D0
!                      !RC( 6)
    xx(3) = 11.0D0
!                      !RC(11)
    yy(1) = zi(nzero) - zi(1)
    yy(2) = zi(nzero) - zi(6)
    yy(3) = zi(nzero) - zi(11)
    bhf(k1,k2) = cautog*ame(k1,k2)*ylag(0.0D0,xx,yy,0,2,3)
  END DO
END DO

END SUBROUTINE corehff
