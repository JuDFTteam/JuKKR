C*==ewald2d.f    processed by SPAG 6.05Rc at 16:01 on 17 May 2004
      SUBROUTINE EWALD2D(LPOT,ALAT,VEC1,VEC2,IQ1,IQ2,RM2,NRMAX,NSHLR,NSR
     &                   ,GN2,NGMAX,NSHLG,NSG,SUM2D,VOL,LASSLD,LMXSPD)
C **********************************************************************
C *                                                                    *
C *   calculation of lattice sums for l .le. 2*lpot :                  *
C *                                                                    *
C *                    ylm( q(i) - q(j) + rm )                         *
C *         sum      ===========================                       *
C *                  | q(i) - q(j) + rm |**(l+1)                       *
C *                                                                    *
C *          - summed over all 2D lattice vectors rm  -                *
C *                                                                    *
C *   ylm       : real spherical harmic to given l,m                   *
C *                                                                    *
C *   The sum is done different in the plane (qi-qj)z = 0              *
C *   and out of the plane. In plane an Ewald procedure similar        *
C *   to the 3d is used and we perform 2 sums (real and reciprocal)    *
C *   the l= 2,4 m=0 terms are calculated with a different method      *
C *                                                                    *
C *   The l=0 term is calculated with a extra factror sqrt(4*pi) this  *
C *   is for transparency reasons (so that the correction terms        *
C *   *r=0,g=0* can be followed in the program)                        *
C *   Literature : lm = (0,0), (1,0) terms :PRB 40, 12164 (1989)       *
C *                                         PRB 49, 2721 (1994)        *
C *                                         PRB 47, 16525 (1993)       *
C *                                         Ziman , p.39-40            *
C *       l=2,4 (m=0) terms are done with recursive diferentiation     *
C *                                                  v. 16.8.99        *
C *       The l multipoles are treated using the expansion             *
C *       for a complex plane wave.                                    *
C *       eq.(33) , M. Weinert, J. Math Phys. 22, 2439 (1981)          *
C *                                                                    *
C *   Final version : 11.01.2000   (No direct sum needed everything    *
C *                                 is done with Ewald method)         *
C *   Programmed by N. Papanikolaou                                    *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Scalar arguments ..
      INTEGER LASSLD,LMXSPD
C      parameter (lassld replaces old l2maxd=2*LPOTD=4*lmaxd)
      INTEGER LPOT,NRMAX,NSHLR,NGMAX,NSHLG,IQ1,IQ2
      DOUBLE PRECISION ALAT,VOL
C     ..
C     .. Array arguments ..
      DOUBLE PRECISION VEC1(3),VEC2(3)
      DOUBLE PRECISION RM2(2,*),GN2(2,*),SUM2D(LMXSPD)
      INTEGER NSR(*),NSG(*)
C     ..
C     .. Local scalars ..
      DOUBLE PRECISION ALPHA,BETA,BOUND,CON,CON1,
     &                 DOT1,DQ1,DQ2,DQ3,DQDOTG,EXPBSQ,FPI,
     &                 G1,G2,G3,GA,LAMDA,
     &                 PI,R,R0,R1,R2,R3,RFAC,S,SIGNRZ,
     &                 STEST0,TPI
      DOUBLE COMPLEX APREFMM,APREFPP,BFAC,CI,FACTEXP,SIMAG
      DOUBLE PRECISION DERFC,EXPONENT,CRIT
      INTEGER I,IM,IR,L,LM,LMAX,LMMAX,M,ICALL,NGMAX1
C     ..
C     .. Local arrays ..
      DOUBLE PRECISION DFAC(0:2*LASSLD+1),G(0:LASSLD),
     &                 GAL(0:LASSLD),GI(0:4),GR(0:4),
     &                 PREF0(0:LASSLD),SIGNRZL(0:LASSLD),YLM(LMXSPD),
     &                 YLMPREF(0:LASSLD),YLMPREF1(0:LASSLD,0:LASSLD)
      DOUBLE COMPLEX CIM(0:LASSLD),EXPONL(0:LASSLD),
     &               PREF2(LASSLD),S0(LMXSPD),
     &               STEST(LMXSPD),STESTNEW(LMXSPD)
C     ..
C     .. External subroutines ..
      EXTERNAL FPLANEG,FPLANER
C     ..
C     .. Intrinsic functions ..
      INTRINSIC DABS,ATAN,DBLE,EXP,SQRT
C     ..
C     .. Data statements
      DATA ICALL /0/
C     ..
C     .. Save statements
      SAVE ICALL,CI,BOUND,PI,FPI,TPI
C     ..................................................................
C
      ICALL = 1
      !ICALL = ICALL + 1
      IF (ICALL.EQ.1) THEN
         CI = (0.0D0,1.0D0)
         BOUND = 1.0D-9
         PI = 4.0D0*ATAN(1.0D0)
         FPI = 4.0D0*PI
         TPI = 2.0D0*PI
      END IF
C
C Factorial
C
      DFAC(0) = 1
      DO L = 1,2*LASSLD + 1
         DFAC(L) = DFAC(L-1)*L
      END DO
      DO L = 0,LASSLD
         PREF0(L) = 0D0
      END DO

      LMAX = 2*LPOT
      LMMAX = (LMAX+1)**2
C
      PREF0(2) = SQRT(5D0/PI)/2D0/2D0
      PREF0(4) = 3D0*SQRT(9D0/PI)/16D0/9D0
C
C --> choose proper splitting parameter
C
      LAMDA = SQRT(PI)/ALAT
C
      DQ1 = (VEC2(1)-VEC1(1))*ALAT   ! SCALE WITH ALAT
      DQ2 = (VEC2(2)-VEC1(2))*ALAT
      DQ3 = (VEC2(3)-VEC1(3))*ALAT
C
C --> Initialise
C
      DO LM = 1,LMMAX
         STEST(LM) = 0D0
         STESTNEW(LM) = 0D0
      END DO
C
C---> Add correction if rz = 0
C
      IF ( DABS(DQ3).LT.1D-6 ) THEN
         STEST(1) = STEST(1) - 2D0*LAMDA/SQRT(PI) - 2D0*SQRT(PI)
     &              /LAMDA/VOL
         STESTNEW(1) = STESTNEW(1) - 2D0*LAMDA/SQRT(PI) - 2D0*SQRT(PI)
     &                 /LAMDA/VOL
C
         IF ( (DQ1*DQ1+DQ2*DQ2).GT.1D-6 ) THEN
            STEST(1) = STEST(1) + 2D0*LAMDA/SQRT(PI)
            STESTNEW(1) = STESTNEW(1) + 2D0*LAMDA/SQRT(PI)
         END IF
      ELSE
C
C---> Add correction if rz<> 0
C
         STEST(1) = STEST(1) - DABS(DQ3)*FPI/2D0/VOL
         STEST(3) = STEST(3) - DABS(DQ3)/DQ3*SQRT(3D0*FPI)/2D0/VOL ! -d/dz
C the correction for higher l vanishes...
      END IF
C
C **********************************************************************
C ******************    I N-P L A N E      M = 0  **********************
C **********************************************************************
      IF ( DABS(DQ3).LT.1D-6 ) THEN
C
C --> Real space sum
C
C ======================================================================
         DO IR = 1,NRMAX
            R1 = DQ1 - RM2(1,IR)
            R2 = DQ2 - RM2(2,IR)
            R3 = 0D0
            R = SQRT(R1*R1+R2*R2)
C-----------------------------------------------------------------------
            IF ( R.GT.1D-8 ) THEN
               ALPHA = LAMDA*R
               CALL FPLANER(ALPHA,GR,R)
               DO L = 0,4
                  LM = L*(L+1) + 1 ! m =0
                  STEST(LM) = STEST(LM) + GR(L)
               END DO
               CALL YMY(R1,R2,R3,R0,YLM,LASSLD)
               CALL GAMFC(ALPHA,G,LASSLD,R)
               YLM(1) = 1D0     ! just definition matter
C
               DO L = 0,LMAX
                  RFAC = G(L)/SQRT(PI)
                  DO M = -L,L
                     LM = L*(L+1) + M + 1
                     STESTNEW(LM) = STESTNEW(LM) + YLM(LM)*RFAC
                  END DO
               END DO
C
               IF ( IR.EQ.(NRMAX-NSR(NSHLR)) ) THEN
C     keep the value before the last shell to test convergence
                  DO L = 0,LMAX
                     DO M = -L,L
                        LM = L*(L+1) + M + 1
                        S0(LM) = STESTNEW(LM)
                     END DO
                  END DO
               END IF
            END IF  ! r <> 0
C-----------------------------------------------------------------------
         END DO                 ! ir loop
C ======================================================================
C
C --> Check convergence
C
         S = 0D0
C ======================================================================
         DO L = 0,LMAX
            DO M = -L,L
               LM = L*(L+1) + M + 1
               STEST0 = ABS(S0(LM)-STESTNEW(LM))
               IF ( S.LT.STEST0 ) S = STEST0
            END DO
         END DO
C ======================================================================
         IF ( S.GT.BOUND ) WRITE (1337,FMT=99001) ABS(S),BOUND,IQ1,IQ2
C
C --> Sum in reciprocal lattice
C
         CON = FPI/2D0/VOL
         ! Find an upper cutoff for G-vectors
         I = 1
         EXPONENT = 1.D0
         CRIT   = 8.D0
         DO WHILE (EXPONENT.LT.CRIT.AND.I.LT.NGMAX)
            I = I + 1
            G1 = GN2(1,I)            ! G vectors are assumed to be sorted with increasing length
            G2 = GN2(2,I)
            GA = SQRT(G1*G1+G2*G2)
            EXPONENT = GA/2.D0/LAMDA  ! If EXPONENT>8., then ERFC(EXPONENT) and EXP(-EXPONENT**2) 
                                      ! (see sub. fplaneg) are below is 1E-27, considered negligible.
         ENDDO
         NGMAX1 = I
!          WRITE(99,FMT='(A7,I8,2E10.2)') 'NGMAX1:',NGMAX1,
!      &                EXP(-EXPONENT**2),ERFC(EXPONENT)

         IF (NGMAX1.GT.NGMAX) STOP 'ewald2d: 1: NGMAX1.GT.NGMAX' ! should never occur
C ======================================================================
         DO IM = 1,NGMAX1
            G1 = GN2(1,IM)
            G2 = GN2(2,IM)
            G3 = 0D0
            GA = SQRT(G1*G1+G2*G2)
C ----------------------------------------------------------------------
            DOT1 = DQ1*G1 + DQ2*G2
            CALL FPLANEG(LAMDA,GI,PREF0,LASSLD,GA,VOL)
            SIMAG = EXP(CI*DOT1)
            DO L = 0,4
               LM = L*(L+1) + 1
               STEST(LM) = STEST(LM) + GI(L)*SIMAG
            END DO
C ----------------------------------------------------------------------
            IF ( GA.GT.1D-6 ) THEN
               CALL YMY(G1,G2,G3,GA,YLM,LASSLD)
               BETA = GA/LAMDA
               EXPBSQ = DERFC(BETA/2D0)
C
               BFAC = CON*SIMAG*EXPBSQ
               STESTNEW(1) = STESTNEW(1) + BFAC/GA
C
               DO L = 0,LMAX
                  IF ( L.NE.0 ) THEN
                     DO M = -L,L
                        LM = L*(L+1) + M + 1
                        STESTNEW(LM) = STESTNEW(LM) + YLM(LM)
     &                                 *BFAC*GA**(L-1)
                     END DO
                  END IF
                  BFAC = BFAC/CI/DBLE(2*L+1)
               END DO
            END IF
C ----------------------------------------------------------------------
            IF ( IM.EQ.(NGMAX1-NSG(NSHLG)) ) THEN
C     keep the value before the last shell to test convergence
               DO LM = 1,LMMAX
                  S0(LM) = STESTNEW(LM)
               END DO
            END IF
         END DO 
C ======================================================================
C
C --> Check convergence
C
         DO LM = 1,LMMAX
            STEST0 = ABS(S0(LM)-STESTNEW(LM))
            IF ( S.LT.STEST0 ) S = STEST0
         END DO
C
C --> Correction due to r=0 term only for DRn = 0
C
C ----------------------------------------------------------------------
         IF ( (DQ1*DQ1+DQ2*DQ2).LT.1D-6 ) THEN
            DO L = 2,4,2
               PREF0(L) = PREF0(L)*DFAC(L)/DFAC(L/2)/(L+1)
            END DO
C
            I = 1
            DO L = 2,4,2
               I = I + 1
               LM = L*(L+1) + 1
               STEST(LM) = STEST(LM) + (-1)**I*2D0/SQRT(PI)*PREF0(L)
     &                     *LAMDA**(L+1)
            END DO
         END IF         
C ----------------------------------------------------------------------
         DO L = 2,4,2
            LM = L*(L+1) + 1
            STESTNEW(LM) = STEST(LM)
         END DO
C --> end of correction
C ----------------------------------------------------------------------
C
         IF ( S.GT.BOUND.AND.NGMAX1.EQ.NGMAX ) 
     &                     WRITE (1337,FMT=99002) ABS(S),BOUND,IQ1,IQ2
C
         DO LM = 1,LMMAX
            STEST(LM) = STESTNEW(LM)
         END DO
C **********************************************************************
      ELSE                      
C **********************************************************************
C ************************* OUT OF THE PLANE ***************************
C **********************************************************************
C
C --> Prepare arrays for speed up
C
         SIGNRZ = DQ3/DABS(DQ3)
         CON1 = TPI/VOL
         DO L = 0,LMAX
            YLMPREF(L) = SQRT((2D0*L+1D0)/FPI)/DFAC(L)*CON1
            SIGNRZL(L) = (-SIGNRZ)**L
            CIM(L) = (-CI)**L
            DO M = 1,L
               YLMPREF1(L,M) = CON1*SQRT(DBLE(2*L+1)/2D0/FPI/DFAC(L+M)
     &                         /DFAC(L-M))
            END DO
         END DO
         YLMPREF(0) = 1D0*CON1
C
C --> Sum in reciprocal space
C
         ! Find an upper cutoff for G-vectors
         I = 1
         EXPONENT = 1.D0
         CRIT   = 60.D0
         DO WHILE (EXPONENT.LT.CRIT.AND.I.LT.NGMAX)
            I = I + 1
            G1 = GN2(1,I)            ! G vectors are assumed to be sorted with increasing length
            G2 = GN2(2,I)
            GA = SQRT(G1*G1+G2*G2)
            EXPONENT = GA*DABS(DQ3)  ! If EXPONENT>60., then EXPBSQ below is 8.7E-27, considered negligible.
         ENDDO
         NGMAX1 = I
         IF (NGMAX1.GT.NGMAX) STOP 'ewald2d: 2: NGMAX1.GT.NGMAX' ! should never occur
!          WRITE(99,FMT='(A7,I8,2E10.2)') 'NGMAX1:',NGMAX1,EXP(-EXPONENT)
C ======================================================================
         DO I = 2,NGMAX1
C
C   Exclude the origin all terms vanish for g=0
C   except the l = 0 components which are treated
C   sepparately. (look at the begining of the sub)
C
            G1 = GN2(1,I)
            G2 = GN2(2,I)
            GA = SQRT(G1*G1+G2*G2)
            DO L = 0,LMAX
               GAL(L) = GA**L
            END DO
C
            EXPBSQ = EXP(-GA*DABS(DQ3))
            DQDOTG = (DQ1*G1+DQ2*G2)
            FACTEXP = EXP(CI*DQDOTG)*EXPBSQ/GA
            EXPONL(0) = 1D0
            DO L = 1,LMAX
               EXPONL(L) = (G1+CI*G2)/GA*EXPONL(L-1)    ! exp(i*m*fi)
            END DO
C
C     In case rz < 0 then multiply by (-1)**(l-m)
C     (M. Weinert J. Math Phys. 22, 2439 (1981) formula 33
C      compare also formula A9)
C
            DO L = 0,LMAX
C     m = 0
               LM = L*(L+1) + 1
               STEST(LM) = STEST(LM) + YLMPREF(L)*GAL(L)*SIGNRZL(L)
     &                     *FACTEXP
C     m <> 0
               DO M = 1,L
                  PREF2(M) = CIM(M)*YLMPREF1(L,M)*GAL(L)
     &                       *SIGNRZ**M*SIGNRZL(L)
C
C Go from the <usual> Jackson Ylm to the ones we use
C
                  APREFPP = (1D0/EXPONL(M)+EXPONL(M))
                  APREFMM = (1D0/EXPONL(M)-EXPONL(M))*CI
C     m > 0
                  LM = L*(L+1) + M + 1
                  STEST(LM) = STEST(LM) + PREF2(M)*APREFPP*FACTEXP
C     m < 0
                  LM = L*(L+1) - M + 1
                  STEST(LM) = STEST(LM) + PREF2(M)*APREFMM*FACTEXP
               END DO
            END DO
C ----------------------------------------------------------------------
            IF ( I.EQ.(NGMAX1-NSG(NSHLG)) ) THEN
C     keep the value before the last shell to test convergence
               DO LM = 1,LMMAX
                  S0(LM) = STEST(LM)
               END DO
            END IF
C ----------------------------------------------------------------------
         END DO
C ======================================================================
C
C --> Check convergence
C
         S = 0D0
         DO LM = 2,LMMAX
            STEST0 = ABS(S0(LM)-STEST(LM))
            IF ( S.LT.STEST0 ) S = STEST0
         END DO
         IF ( S.GT.BOUND.AND.NGMAX1.EQ.NGMAX ) 
     &                     WRITE (1337,FMT=99003) ABS(S),BOUND,IQ1,IQ2
      END IF 
C **********************************************************************
C
      DO LM = 1,LMMAX
         IF ( ABS(DIMAG(STEST(LM))).GT.BOUND ) THEN
            WRITE (6,*) ' ERROR: Imaginary contribution',
     &                  ' to REAL lattice sum',DIMAG(STEST(LM)),BOUND
            STOP
         END IF
         SUM2D(LM) = DBLE(STEST(LM))
         STEST(LM) = 0.0D0
      END DO
C
      SUM2D(1) = SUM2D(1)/SQRT(FPI)
C
99001 FORMAT (5X,
     &     'WARNING 1 : Convergence of 2D-sum is ',1P,E9.2,' > ',E9.2,
     &     'LAYER PAIR',2I6,/,15X,
     &     'You should use more lattice vectors (RMAX)')
99002 FORMAT (5X,
     &     'WARNING 2 : Convergence of 2D-sum is ',1P,E9.2,' > ',E9.2,
     &     'LAYER PAIR',2I6,/,15X,
     &     'You should use more lattice vectors (GMAX)')
99003 FORMAT (5X,
     &     'WARNING 3 : Convergence of 2D-sum is ',1P,E9.2,' > ',E9.2,
     &     'LAYER PAIR',2I6,/,15X,
     &     'You should use more lattice vectors (GMAX)')
      END
