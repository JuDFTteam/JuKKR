! **********************************************************************
      SUBROUTINE GFREE13(RDIFF,E0,GMLL,DGMLL,CLEB,ICLEB,LOFLM,IEND)
! **********************************************************************
!     .. Parameters ..
      IMPLICIT NONE
      include 'inc.p'
      INTEGER LMAX
      PARAMETER (LMAX=LMAXD)
      INTEGER LMGF0D
      PARAMETER (LMGF0D= (LMAX+1)**2)
      INTEGER LMAX2P,LMX2SQ
      PARAMETER (LMAX2P=LMAX*2+1,LMX2SQ=LMAX2P**2)
      DOUBLE COMPLEX CZERO,CI
      PARAMETER (CZERO= (0.0D0,0.0D0),CI= (0.0D0,1.0D0))
!     ..
!     .. Scalar Arguments ..
      DOUBLE COMPLEX E0
      INTEGER IEND
!     ..
!     .. Array Arguments ..
      DOUBLE COMPLEX GMLL(LMGF0D,LMGF0D)
      DOUBLE COMPLEX DGMLL(LMGF0D,LMGF0D) ! LLY
      DOUBLE PRECISION CLEB(NCLEB),RDIFF(*)
      INTEGER ICLEB(NCLEB,4),LOFLM(*)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION FPI,PI,RABS,RFPI,X,Y,Z
      INTEGER IFAC,J,LM1,LM2,LM3
!     ..
!     .. Local Arrays ..
      DOUBLE COMPLEX BL(LMAX2P),HL(LMAX2P),HYL(LMX2SQ),NL(LMAX2P)
      DOUBLE COMPLEX DHL(LMAX2P),DHYL(LMX2SQ) ! LLY
      DOUBLE PRECISION YL(LMX2SQ)
      INTEGER LF(LMX2SQ)
!     ..
!     .. External Subroutines ..
      EXTERNAL BESHAN,YMY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
!     ..
      PI = 4.D0*ATAN(1.D0)
      FPI = 4.D0*PI
      RFPI = SQRT(FPI)
!-----------------------------------------------------------------------
!---- CALCULATION OF FREE ELECTRON GREEN'S FUNCTION :  G(M)LL'(E0)
!-----------------------------------------------------------------------
!     Also:
!---- derivative of free electron green function matrix elements
!     returned in array DGMLL=dGMLL / dE
!
!     the analytical formula for the derivative of spherical Hankel
!     functions is used:
!
!     d                     l+1
!     --  h (x) = h   (x) - --- h (x)   
!     dx   l       l-1       x   l
!
!     which for x = sqrt(E0)*r leads to
!
!      d                       r           rl
!     --- ( sqrt(E0) h (x) ) = - h   (x) - -- h (x) )
!     dE0             l        2  l-1      2x  l
!
!     Ported from KKRnano by Phivos Mavropoulos 10.10.2013
!-----------------------------------------------------------------------
      DO 10 LM1 = 1,LMX2SQ
        LF(LM1) = LOFLM(LM1) + 1
   10 CONTINUE
      X = RDIFF(1)
      Y = RDIFF(2)
      Z = RDIFF(3)
      CALL YMY(X,Y,Z,RABS,YL,LMAX*2)
      CALL BESHAN(HL,BL,NL,SQRT(E0)*RABS,LMAX*2)

      ! Derivatives of Hankel functions ! LLY
      DHL(1) = 0.5D0*CI*RABS*HL(1)
      DO LM1 = 2,LMAX2P
         DHL(LM1) = 0.5D0*(RABS*HL(LM1-1)-(LM1-1)*HL(LM1)/SQRT(E0))
      END DO 

      DO 20 LM1 = 1,LMX2SQ
        HYL(LM1) = -FPI*CI*SQRT(E0)*YL(LM1)*HL(LF(LM1))
        DHYL(LM1) = -FPI*CI*YL(LM1)*DHL(LF(LM1))   ! LLY
   20 CONTINUE


      DO 40 LM1 = 1,LMGF0D
        GMLL(LM1,LM1) = HYL(1)/RFPI
        DGMLL(LM1,LM1) = DHYL(1)/RFPI  ! LLY
        DO 30 LM2 = 1,LM1 - 1
          GMLL(LM1,LM2) = CZERO
          DGMLL(LM1,LM2) = CZERO    ! LLY
   30   CONTINUE
   40 CONTINUE

      DO 50 J = 1,IEND
        LM1 = ICLEB(J,1)
        LM2 = ICLEB(J,2)
        LM3 = ICLEB(J,3)
        GMLL(LM1,LM2) = GMLL(LM1,LM2) + CLEB(J)*HYL(LM3)
        DGMLL(LM1,LM2) = DGMLL(LM1,LM2) + CLEB(J)*DHYL(LM3)  ! LLY
   50 CONTINUE
      DO 70 LM1 = 1,LMGF0D
        DO 60 LM2 = 1,LM1 - 1
          IFAC = (-1)** (LOFLM(LM1)+LOFLM(LM2))
          GMLL(LM2,LM1) = IFAC*GMLL(LM1,LM2)
          DGMLL(LM2,LM1) = IFAC*DGMLL(LM1,LM2)  ! LLY
   60   CONTINUE
   70 CONTINUE
      RETURN

      END
