      MODULE MOD_JMTRX

      CONTAINS

      SUBROUTINE JMTRX(X1,X2,X3,E,LMAX,JMAT,LCALL)
      use mod_gauntharmonics, only: gauntcoeff
      use mod_ymy
      use mod_bessel1
      implicit none
! C ****************************************************************
! C * This subroutine calculates the transformation matrix Ull'(E) *
! C * for the Green's function in the shifted position. According  *
! C * N. Stefanou et al PRB 36 6372 (1987) eq.(6).                 *
! C * Lmax is the maximum l cutoff for the matrix Ull':=JMAT       *
! C * Since the Gaunt coefficients are used, subroutines gaunt and *
! C * Gaunt2 have to be called first to set up the common block    *
! C * GAUNTC.                                                      *
! C ****************************************************************
! C
! C     .. Parameters ..
! c
! c---> attention : ncleb is an empirical factor - it has to be optimized
! c
!       include 'parameters.file'
! c      INTEGER LMX,LPOTD,LMAXD
! c      PARAMETER (LMAXD=3,LPOTD=6,LMX=LMAXD+1)
!       INTEGER LMMAXD,LM2D,LMPOTD
!       PARAMETER (LMMAXD= (LMAXD+1)**2,LM2D= (2*LMAXD+1)**2,
!      +          LMPOTD= (LPOTD+1)**2)
!       INTEGER LMMAX2D
!       PARAMETER (TWOLMAX=2*LMAXD,LMMAX2D=(TWOLMAX+1)**2)
!       INTEGER N
!       PARAMETER (N=4*LMAXD)
!       INTEGER NCLEB
!       PARAMETER (NCLEB=LMPOTD*LMMAXD)
! C     ..
! C     .. Scalar Arguments ..
! C
!       INTEGER LMAX,IA,L1,L2,L3,LM1,LM2,LM3
!       INTEGER LMMAX,I,LM,L
      COMPLEX*16 ARG,E,CLLL1
      REAL*8 SABS,FOURP,X1,X2,X3,C0LL,X123
      LOGICAL  LCALL
! C
! C     .. Array arguments ..
! C
      REAL*8 YLM((2*LMAX+1)**2)
      COMPLEX*16 BJ(0:2*LMAX),H(0:2*LMAX),Y(0:2*LMAX)
      COMPLEX*16 JMAT((LMAX+1)**2,(LMAX+1)**2)
! C
! C     .. Intrinsic functions ..
! C
      INTRINSIC DSQRT,DATAN,CDSQRT
! C
! C     .. External subroutines ..
! C
!       EXTERNAL YMY,BESSEL1
! C
! C     .. Scalars in Common ..
! C
! !       INTEGER IEND
! C     ..
! C     .. Arrays in Common ..
! C
! !       REAL*8 CLEB(NCLEB,2)
! !       INTEGER ICLEB(NCLEB,4),JEND(LMPOTD,0:LMAXD,0:LMAXD),LOFLM(LM2D)
! C
!       COMMON /GAUNTC/CLEB,LOFLM,ICLEB,IEND,JEND
      DOUBLE COMPLEX,parameter ::  EI=   (0.0D0,1.0D0)
      DOUBLE COMPLEX,parameter ::  CZERO=(0.0D0,0.0D0)
      INTEGER :: LMMAX,LMAX,LM,LM1,IA,LM2,LM3
      INTEGER :: l1,l2,l3

      FOURP=16.D0*DATAN(1.D0)
      C0LL=1.D0/SQRT(FOURP)
      LMMAX=(LMAX+1)**2

      X123=SQRT(X1*X1+X2*X2+X3*X3)

    IF (X123>1.D-10) THEN

      CALL YMY(X1,X2,X3,SABS,YLM,2*LMAX)
      ARG=SQRT(E)*SABS
      CALL BESSEL1(BJ,Y,H,ARG,2*LMAX,2*LMAX,.TRUE.,.TRUE.,.TRUE.,LCALL)
      DO LM=1,LMMAX
        DO LM1=1,LMMAX
          JMAT(LM,LM1)=CZERO
        END DO
      END DO

      DO IA=1,gauntcoeff(LMAX)%IEND
        LM1=gauntcoeff(LMAX)%ICLEB(IA,1)
        LM2=gauntcoeff(LMAX)%ICLEB(IA,2)
        LM3=gauntcoeff(LMAX)%ICLEB(IA,3)
        L1=gauntcoeff(LMAX)%LOFLM(LM1)
        L2=gauntcoeff(LMAX)%LOFLM(LM2)
        L3=gauntcoeff(LMAX)%LOFLM(LM3)
!       calculate jmat for lm1.ge.lm2
          write(*,'(3I,10F)') lm1,lm2,lm3,gauntcoeff(LMAX)%CLEB(IA,1),bj(l3),ylm(lm3)
        CLLL1=FOURP*gauntcoeff(LMAX)%CLEB(IA,1)*(EI)**(L1+L3-L2)
        JMAT(LM1,LM2)=JMAT(LM1,LM2)+CLLL1*BJ(L3)*YLM(LM3)
      END DO
!  add the l=0 component in the diagonal elements.
      DO LM1=1,LMMAX
        JMAT(LM1,LM1)=JMAT(LM1,LM1)+FOURP*C0LL*BJ(0)*YLM(1)
      END DO
!   Create the rest of the matrix by symmetry relation.
      DO LM1=1,LMMAX
        DO LM2=1,LM1-1
          L1=gauntcoeff(LMAX)%LOFLM(LM1)
          L2=gauntcoeff(LMAX)%LOFLM(LM2)
          JMAT(LM2,LM1)=JMAT(LM1,LM2)*(-1.D0)**(L1+L2)
        END DO
      END DO

    ELSE ! r123<10e-10

      DO LM=1,LMMAX
        DO LM1=1,LMMAX
          JMAT(LM,LM1)=CZERO
          JMAT(LM,LM)=CMPLX(1.0D0,0.0D0)
        END DO
      END DO

    END IF

    END SUBROUTINE
    END MODULE MOD_JMTRX
