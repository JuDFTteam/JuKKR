! OpenMP parallelised, needs threadsafe erfcex, gamfc and ymy E.R.
SUBROUTINE STRMAT(ALAT,LPOT,NAEZ,NGMAX,NRMAX,NSG,NSR,NSHLG,NSHLR, &
GN,RM,QI0,SMAT,VOL,LASSLD,LMXSPD,NAEZD,I1)
  ! **********************************************************************
  ! *                                                                    *
  ! *  calculation of lattice sums for l .le. 2*lpot :                   *
  ! *                                                                    *
  ! *                   ylm( q(i) - q(j) - rm )                          *
  ! *        sum      ===========================                        *
  ! *                 | q(i) - q(j) - rm |**(l+1)                        *
  ! *                                                                    *
  ! *         - summed over all lattice vectors rm  -                    *
  ! *                                                                    *
  ! *  ylm       : real spherical harmic to given l,m                    *
  ! *  q(i),q(j) : basis vectors of the unit cell                        *
  ! *                                                                    *
  ! *  in the case of i = j, rm = 0 is omitted.                          *
  ! *                                                                    *
  ! *  the ewald method is used to perform the lattice summations        *
  ! *  the splitting parameter lamda is set equal sqrt(pi)/alat          *
  ! *  (alat is the lattice constant) .                                  *
  ! *                                                                    *
  ! *  if the contribution of the last shell of the direct and the       *
  ! *  reciprocal lattice is greater than 1.0e-8 a message is written    *
  ! *                                                                    *
  ! *                                    b.drittler may 1989             *
  ! *                                                                    *
  ! *  Dimension of arrays GN,RM changed from (4,*) to (3,*), the 4th    *
  ! *  one not being used (see also lattice3d)     v.popescu May 2004    *
  ! *                                                                    *
  ! **********************************************************************
  IMPLICIT NONE
  !     ..
  !     .. Parameters ..
  DOUBLE COMPLEX CI
  PARAMETER ( CI = (0D0,1D0) )
  DOUBLE PRECISION BOUND
  PARAMETER ( BOUND=1D-8 )
  !     ..
  !     .. Scalar arguments ..
  DOUBLE PRECISION ALAT,VOL
  INTEGER LPOT,NAEZ,NGMAX,NRMAX,NSHLG,NSHLR
  INTEGER LASSLD,LMXSPD,NAEZD
  !     ..
  !     .. Array arguments ..
  DOUBLE PRECISION  GN(3,*),QI0(3,*),RM(3,*),SMAT(LMXSPD,*)
  INTEGER NSG(*),NSR(*)
  !     ..
  !     .. Local scalars ..
  DOUBLE COMPLEX BFAC
  DOUBLE PRECISION ALPHA,BETA,DQ1,DQ2,DQ3,DQDOTG,EXPBSQ,FPI, &
  G1,G2,G3,GA,LAMDA,PI,R,R1,R2,R3,RFAC,S
  DOUBLE PRECISION DBLE
  INTEGER I,I1,I2,IT,L,LM,LMX,LMXSP,M,NGE,NGS,NRE,NRS,NSTART
  !     ..
  !     .. Local arrays ..
  DOUBLE COMPLEX STEST(LMXSPD)
  DOUBLE PRECISION G(0:LASSLD),YLM(LMXSPD) ! QI(3,NAEZD) get rid of this array
  !     ..
  !     .. External subroutines ..
  EXTERNAL GAMFC,YMY
  !     ..
  !     .. Intrinsic functions ..
  INTRINSIC ATAN,ABS,DBLE,EXP,SQRT
  !     ..................................................................

  LMX = 2*LPOT
  LMXSP = (LMX+1)*(LMX+1)
  PI = 4.0D0*ATAN(1.0D0)
  FPI = 4.0D0*PI

  ! --> choose proper splitting parameter

  LAMDA = SQRT(PI)/ALAT

  ! --> loop over atoms per unit cell -- scale basis atoms with alat

  !DO I2 = 1,NAEZ
  !   QI(1,I2) = QI0(1,I2)*ALAT
  !   QI(2,I2) = QI0(2,I2)*ALAT
  !   QI(3,I2) = QI0(3,I2)*ALAT
  !END DO

  ! **********************************************************************
  !$omp parallel do private(I2,DQ1,DQ2,DQ3,STEST,LM,NSTART,IT, &
  !$                        NRS,NGS,NRE,NGE,I,R1,R2,R3, &
  !$                        YLM,R,ALPHA,G,RFAC,L,M, &
  !$                        G1,G2,G3,GA,BETA,EXPBSQ,DQDOTG,BFAC,S)
  DO I2 = 1,NAEZ
  !======================================================================
    DQ1 = (QI0(1,I1) - QI0(1,I2)) * ALAT
    DQ2 = (QI0(2,I1) - QI0(2,I2)) * ALAT
    DQ3 = (QI0(3,I1) - QI0(3,I2)) * ALAT

    STEST(1) = -SQRT(FPI)/VOL/(4D0*LAMDA*LAMDA)
    DO LM = 2,LMXSP
      STEST(LM) = 0.0D0
    END DO

    ! --> exclude the origine and add correction if i1.eq.i2

    IF ( I1.EQ.I2 ) THEN
      STEST(1) = STEST(1) - LAMDA/PI
      NSTART = 2
    ELSE
      NSTART = 1
    END IF
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! --> loop first over n-1 shells of real and reciprocal lattice - then
    !     add the contribution of the last shells to see convergence

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DO IT = 1,2
      IF ( IT.EQ.1 ) THEN
        NRS = NSTART
        NGS = 2
        NRE = NRMAX - NSR(NSHLR)
        NGE = NGMAX - NSG(NSHLG)
      ELSE
        NRS = NRE + 1
        NGS = NGE + 1
        NRE = NRMAX
        NGE = NGMAX
      END IF

      ! --> sum over real lattice

      ! ---------------------------------------------------------------------
      DO I = NRS,NRE
        R1 = DQ1 - RM(1,I)
        R2 = DQ2 - RM(2,I)
        R3 = DQ3 - RM(3,I)

        CALL YMY(R1,R2,R3,R,YLM,LMX)
        ALPHA = LAMDA*R
        CALL GAMFC(ALPHA,G,LMX,R)

        DO L = 0,LMX
          RFAC = G(L)/SQRT(PI)
          DO M = -L,L
            LM = L*(L+1) + M + 1
            STEST(LM) = STEST(LM) + YLM(LM)*RFAC
          END DO
        END DO

      END DO
      ! ---------------------------------------------------------------------

      ! --> sum over reciprocal lattice

      ! ---------------------------------------------------------------------
      DO I = NGS,NGE
        G1 = GN(1,I)
        G2 = GN(2,I)
        G3 = GN(3,I)

        CALL YMY(G1,G2,G3,GA,YLM,LMX)
        BETA = GA/LAMDA
        EXPBSQ = EXP(BETA*BETA/4.0D0)
        DQDOTG = DQ1*G1 + DQ2*G2 + DQ3*G3

        BFAC = FPI*EXP(CI*DQDOTG)/(GA*GA*EXPBSQ*VOL)

        DO L = 0,LMX
          DO M = -L,L
            LM = L*(L+1) + M + 1
            STEST(LM) = STEST(LM) + YLM(LM)*BFAC
          END DO
          BFAC = BFAC*GA/CI/DBLE(2*L+1)
        END DO
      END DO
      ! ---------------------------------------------------------------------
      IF ( IT.EQ.1 ) THEN
        DO LM = 1,LMXSP
          IF ( ABS(DIMAG(STEST(LM))).GT.BOUND ) THEN
            WRITE (6,*) ' ERROR: Imaginary contribution', &
            ' to REAL lattice sum'
            STOP
          END IF
          SMAT(LM,I2) = DBLE(STEST(LM))
          STEST(LM) = 0.0D0
        END DO
      ELSE

        ! --> test convergence

        DO LM = 1,LMXSP
          S = DBLE(STEST(LM))
          SMAT(LM,I2) = SMAT(LM,I2) + S
          !IF (2.LT.1.AND. ABS(S).GT.BOUND ) WRITE (6,FMT=99001) I1,I2, &
          !LM,ABS(S)
        END DO
      END IF
    ! ---------------------------------------------------------------------
    END DO
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  END DO   ! I2
  !$omp end parallel do
  ! **********************************************************************

99001 FORMAT (5X,'WARNING : Convergence of SMAT(',I2,',',I2,') ', &
  ' for LMXSP =',I3,' is ',1P,D8.2,' > 1D-8',/,15X, &
  'You should use more lattice vectors (RMAX/GMAX)')
END
