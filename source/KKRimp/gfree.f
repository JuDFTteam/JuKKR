!------------------------------------------------------------------------------------
!> Summary: This module is intended for the calculation of the free electron Green function
!> Author: People who wrote it
!> A More detailed explanation with the math, concepts, etc necessary to understand the routine
!>
!> One can write Latex comments like this \(i\hbar\frac{\partial \psi}{\partial t}=-\mathcal{H}\psi\)
!> or add labeled equations using the standard latex way
!> \begin{equation}
!> \mathbf{A} \mathbf{v}= \eta\mathbf{v}
!> \end{equation}
!> **FORd** also accepts markdown style so you can _write with style_ 
!> 
!> **IMPORTANT**
!> The JM-KKR follows the coding conventions noted in this example, one that is
!> not obvious is that **each level of indentation consists of two spaces**. Please keep this:
!> _do it for the children_.
!> So please keep the conventions.
! These are special boxes for ford, notice that this comment does not appear in the html file.
! These boxes contain important information and should be added when necessary. ALWAYS remember to close the box
! BEFORE oppening a new one or they will be nested.
!------------------------------------------------------------------------------------

      MODULE mod_gfree
      CONTAINS

!-------------------------------------------------------------------------------
!> Summary: Calculate the free electron Green function
!> Author: Who wrote this subroutine
!> Category: Greenfunction, reference-system
!> Deprecated: TRUE ! This needs to be set to True for deprecated subroutines
!> A More detailed explanation with the math, concepts, etc necessary to understand the routine
!-------------------------------------------------------------------------------


c **********************************************************************
      SUBROUTINE GFREE(RDIFF,E0,GMLL,CLEB,ICLEB,LOFLM,IEND,
     +                 NCLEB,LMAX,LMGF0D,LMAX2P,LMX2SQ)
       USE mod_beshan
       USE mod_ymy
      IMPLICIT NONE
c **********************************************************************
C     .. Parameters ..
!       include 'inc.p'
       INTEGER LMAX
!       PARAMETER (LMAX=LMAXD)
       INTEGER NCLEB
      INTEGER LMGF0D
!       PARAMETER (LMGF0D= (LMAX+1)**2)
      INTEGER LMAX2P,LMX2SQ
!       PARAMETER (LMAX2P=LMAX*2+1,LMX2SQ=LMAX2P**2)
       DOUBLE COMPLEX CZERO,CI
       PARAMETER (CZERO= (0.0D0,0.0D0),CI= (0.0D0,1.0D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX E0
      INTEGER IEND
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX GMLL(LMGF0D,LMGF0D)
      DOUBLE PRECISION CLEB(NCLEB),RDIFF(:)
      INTEGER ICLEB(NCLEB,4),LOFLM(:)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION FPI,PI,RABS,RFPI,X,Y,Z
      INTEGER IFAC,J,LM1,LM2,LM3
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX BL(LMAX2P),HL(LMAX2P),HYL(LMX2SQ),NL(LMAX2P)
      DOUBLE PRECISION YL(LMX2SQ)
      INTEGER LF(LMX2SQ)
C     ..
C     .. External Subroutines ..
!       EXTERNAL YMY
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
C     ..
      PI = 4.D0*ATAN(1.D0)
      FPI = 4.D0*PI
      RFPI = SQRT(FPI)
C-----------------------------------------------------------------------
C---- CALCULATION OF FREE ELECTRON GREEN'S FUNCTION :  G(M)LL'(E0)
C-----------------------------------------------------------------------
      DO 10 LM1 = 1,LMX2SQ
        LF(LM1) = LOFLM(LM1) + 1
   10 CONTINUE
      X = RDIFF(1)
      Y = RDIFF(2)
      Z = RDIFF(3)
      CALL YMY(X,Y,Z,RABS,YL,LMAX*2)
      CALL BESHAN(HL,BL,NL,SQRT(E0)*RABS,LMAX*2)
      DO 20 LM1 = 1,LMX2SQ
        HYL(LM1) = -FPI*CI*SQRT(E0)*YL(LM1)*HL(LF(LM1))
   20 CONTINUE
      DO 40 LM1 = 1,LMGF0D
        GMLL(LM1,LM1) = HYL(1)/RFPI
        DO 30 LM2 = 1,LM1 - 1
          GMLL(LM1,LM2) = CZERO
   30   CONTINUE
   40 CONTINUE
      DO 50 J = 1,IEND
        LM1 = ICLEB(J,1)
        LM2 = ICLEB(J,2)
        LM3 = ICLEB(J,3)
        IF (LM1 <= LMGF0D .AND. LM2 <= LMGF0D ) THEN !David2009
          GMLL(LM1,LM2) = GMLL(LM1,LM2) + CLEB(J)*HYL(LM3)
        END IF
   50 CONTINUE
      DO 70 LM1 = 1,LMGF0D
        DO 60 LM2 = 1,LM1 - 1
          IFAC = (-1)** (LOFLM(LM1)+LOFLM(LM2))
          GMLL(LM2,LM1) = IFAC*GMLL(LM1,LM2)
   60   CONTINUE
   70 CONTINUE
      RETURN

      END SUBROUTINE
      END MODULE mod_GFREE
