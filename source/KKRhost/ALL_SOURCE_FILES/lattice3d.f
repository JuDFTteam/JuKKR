C*==lattice3d.f    processed by SPAG 6.05Rc at 19:02 on 18 May 2004
      SUBROUTINE LATTICE3D(ALAT,BRAVAIS,RECBV,NGMAX,NRMAX,NSHLG,NSHLR,
     +                     NSG,NSR,GN,RM,RMAX,GMAX,IPRINT,NMAXD,ISHLD)
C **********************************************************************
C *                                                                    *
C *  generate lattice vectors of direct and reciprocal space from      *
C *  basic translation vectors br                                      *
C *                                                                    *
C *  alat            : lattice constant                                *
C *  br(i,j)         : i=x,y,z j= 1,2,3 bravais vectors                *
C *                    *** in a.u. ****                                *
C *  rmax            : maximum radius in real space        (input)     *
C *  gmax            : maximum radius in reciprocal space  (input)     *
C *  ngmax           : Number of reciprocal lattice vectors            *
C *  gn(3,nmaxd)     : x,y,z   of reciprocal lattice vectors           *
C *  nrmax           : Number of real lattice vectors                  *
C *  rm(3,nmaxd)     : x,y,z  of real space vectors                    *
C *  nshlg           : shells in reciprocal space                      *
C *  nshlr           : shells in real space                            *
C *  nsg,nsr         : integer arrays, number of atoms in each shell   *
C *                                                                    *
C *  Dimension of arrays GN,RM changed from (4,*) to (3,*), the 4th    *
C *  one it is used only locally (GNR/RMR)       v.popescu May 2004    *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Scalar arguments ..
      INTEGER IPRINT,NGMAX,NRMAX,NSHLG,NSHLR,NMAXD,ISHLD
      DOUBLE PRECISION ALAT,RMAX,GMAX
C     ..
C     .. Array arguments ..
      DOUBLE PRECISION BRAVAIS(3,3),RECBV(3,3)
      DOUBLE PRECISION GN(3,NMAXD),RM(3,NMAXD)
      INTEGER NSG(ISHLD),NSR(ISHLD)
C     ..
C     .. Local scalars ..
      DOUBLE PRECISION A,ABSGM,ABSRM,AG,AR,B,C,DA,DB,GX,GY,GZ,PI,
     &                 RX,RY,RZ,VMIN
      INTEGER I,K,L,M,N,N1,NG,NR,NSH,NSHL,NUMG,NUMGH,NUMR,
     &        NUMRH
      DOUBLE PRECISION DBLE
      INTEGER IDINT
C     ..
C     .. Local arrays ..
      DOUBLE PRECISION ABSG(3),ABSR(3),BG(3,3),BR(3,3),CJ(4,NMAXD)
      DOUBLE PRECISION GNR(NMAXD),RMR(NMAXD)
C     ..
C     .. Intrinsic functions ..
      INTRINSIC ABS,ATAN,DBLE,IDINT,MAX,MOD,SQRT
C     ..
C     .. External subroutines ..
      EXTERNAL IOINPUT
C     ..................................................................
C
      PI = 4.0D0*ATAN(1.0D0)
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (1337,'(5X,2A,/)') '< LATTICE3D > : ',
     &                    'generating direct/reciprocal lattice vectors'
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
      RMAX = RMAX*ALAT
      GMAX = GMAX/ALAT
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (1337,FMT=99001) RMAX,GMAX
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C ======================================================================
C
C --> basic trans. vectors and basis vectors
C
      DO I = 1,3
         BR(1,I) = BRAVAIS(1,I)*ALAT
         BR(2,I) = BRAVAIS(2,I)*ALAT
         BR(3,I) = BRAVAIS(3,I)*ALAT
      END DO
C ======================================================================
C
C --> generate primitive vectors BG of reciprocal space
C
      DO I = 1,3
         BG(1,I) = RECBV(1,I)*2D0*PI/ALAT
         BG(2,I) = RECBV(2,I)*2D0*PI/ALAT
         BG(3,I) = RECBV(3,I)*2D0*PI/ALAT
      END DO
C ======================================================================
C
C --> estimate no. of lattice vectors
C
      DO I = 1,3
         ABSR(I) = SQRT(BR(1,I)**2+BR(2,I)**2+BR(3,I)**2)
         ABSG(I) = SQRT(BG(1,I)**2+BG(2,I)**2+BG(3,I)**2)
      END DO
C
      ABSRM = MAX(ABSR(1),ABSR(2),ABSR(3))
      ABSGM = MAX(ABSG(1),ABSG(2),ABSG(3))
      ABSRM = 2.0D0*PI/ABSRM
      ABSGM = 2.0D0*PI/ABSGM
      NUMR = 2*(IDINT(RMAX/ABSGM)+1) + 1
      NUMG = 2*(IDINT(GMAX/ABSRM)+1) + 1
      NUMRH = NUMR/2 + 1
      NUMGH = NUMG/2 + 1
C
C **********************************************************************
C                 generate lattice vectors of real space
C **********************************************************************
C
      NR = 0
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO L = 1,NUMR
         A = DBLE(L-NUMRH)
         DO M = 1,NUMR
            B = DBLE(M-NUMRH)
            DO N = 1,NUMR
               C = DBLE(N-NUMRH)
C ----------------------------------------------------------------------
               RX = A*BR(1,1) + B*BR(1,2) + C*BR(1,3)
               RY = A*BR(2,1) + B*BR(2,2) + C*BR(2,3)
               RZ = A*BR(3,1) + B*BR(3,2) + C*BR(3,3)
               AR = SQRT(RX*RX+RY*RY+RZ*RZ)
C ----------------------------------------------------------------------
               IF ( AR.LE.RMAX ) THEN
                  NR = NR + 1
                  IF ( NR.GT.NMAXD ) THEN
                     WRITE (6,*) 
     &                      ' ERROR: Dimension NMAXD in inc.p too small'
     &                      ,NR,NMAXD
                     STOP
                  END IF
                  CJ(1,NR) = RX
                  CJ(2,NR) = RY
                  CJ(3,NR) = RZ
                  CJ(4,NR) = AR
               END IF
            END DO
         END DO
      END DO
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      NRMAX = NR
C ======================================================================
C
C --> sort vectors in order of increasing absolute value
C
      DA = 1.D-06
      NSH = 0
      NSHL = -1
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO K = 1,NR
         VMIN = RMAX + 1.0D0
         DO N = 1,NR
            IF ( CJ(4,N)-VMIN.LT.0D0 ) THEN
               VMIN = CJ(4,N)
               N1 = N
            END IF
         END DO
C
         NSHL = NSHL + 1
         RM(1,K) = CJ(1,N1)
         RM(2,K) = CJ(2,N1)
         RM(3,K) = CJ(3,N1)
         RMR(K)  = CJ(4,N1)
         DB = VMIN
C ----------------------------------------------------------------------
         IF ( DB.GT.DA+1.D-06 ) THEN
            NSH = NSH + 1
            IF ( NSH.GT.ISHLD ) THEN
               WRITE (6,*) ' ERROR: Dimension ISHLD in inc.p too small',
     &                     NSH,ISHLD
               STOP
            END IF
C
            NSR(NSH) = NSHL
            NSHL = 0
            DA = DB
         END IF
C ----------------------------------------------------------------------
         CJ(4,N1) = RMAX + 1.0D0
      END DO
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      NSH = NSH + 1
      NSHL = NSHL + 1
      IF ( NSH.GT.ISHLD ) THEN
         WRITE (6,*) ' ERROR: Dimension ISHLD in inc.p too small',NSH,
     &               ISHLD
         STOP
      END IF
C
      NSR(NSH) = NSHL
      NSHLR = NSH
      IF ( NSHLR.LE.1 ) STOP ' ERROR: cut-off radius RMAX too small '
C
C **********************************************************************
C
C **********************************************************************
C                 generate lattice vectors of real space
C **********************************************************************
C
      NG = 0
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO L = 1,NUMG
         A = DBLE(L-NUMGH)
         DO M = 1,NUMG
            B = DBLE(M-NUMGH)
            DO N = 1,NUMG
               C = DBLE(N-NUMGH)
C ----------------------------------------------------------------------
               GX = A*BG(1,1) + B*BG(1,2) + C*BG(1,3)
               GY = A*BG(2,1) + B*BG(2,2) + C*BG(2,3)
               GZ = A*BG(3,1) + B*BG(3,2) + C*BG(3,3)
               AG = SQRT(GX*GX+GY*GY+GZ*GZ)
C ----------------------------------------------------------------------
               IF ( AG.LE.GMAX ) THEN
                  NG = NG + 1
                  IF ( NG.GT.NMAXD ) THEN
                     WRITE (6,*) 
     &                      ' ERROR: Dimension NMAXD in inc.p too small'
     &                      ,NG,NMAXD
                     STOP
                  END IF
                  CJ(1,NG) = GX
                  CJ(2,NG) = GY
                  CJ(3,NG) = GZ
                  CJ(4,NG) = AG
               END IF
            END DO
         END DO
      END DO
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      NGMAX = NG
C ======================================================================
C
C --> sort vectors in order of increasing abs. value
C
      DA = 1.D-06
      NSH = 0
      NSHL = -1
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO K = 1,NG
         VMIN = GMAX + 1.0D0
         DO N = 1,NG
            IF ( CJ(4,N).LT.VMIN ) THEN
               VMIN = CJ(4,N)
               N1 = N
            END IF
         END DO
C
         NSHL = NSHL + 1
         GN(1,K) = CJ(1,N1)
         GN(2,K) = CJ(2,N1)
         GN(3,K) = CJ(3,N1)
         GNR(K)  = CJ(4,N1)
         DB = VMIN
C ----------------------------------------------------------------------
         IF ( DB.GT.DA+1.D-07 ) THEN
            NSH = NSH + 1
            IF ( NSH.GT.ISHLD ) THEN
               WRITE (6,*) ' ERROR: Dimension ISHLD in inc.p too small',
     &                     NSH,ISHLD
               STOP
            END IF
C
            NSG(NSH) = NSHL
            NSHL = 0
            DA = DB
         END IF
C ----------------------------------------------------------------------
         CJ(4,N1) = GMAX + 1.0D0
      END DO
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      NSH = NSH + 1
      NSHL = NSHL + 1
      IF ( NSH.GT.ISHLD ) THEN
         WRITE (6,*) ' ERROR: Dimension ISHLD in inc.p too small',NSH,
     &               ISHLD
         STOP
      END IF
C
      NSG(NSH) = NSHL
      NSHLG = NSH
      IF ( NSHLG.LE.1 ) STOP ' ERROR: cut-off radius GMAX too small '
C **********************************************************************
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (1337,FMT=99002)
      WRITE (1337,FMT=99003) 'Direct  lattice',NRMAX,NSHLR,RMR(NRMAX)
      WRITE (1337,FMT=99003) 'Recipr. lattice',NGMAX,NSHLG,GNR(NGMAX)
      WRITE (1337,FMT=99004)
C
      IF ( IPRINT.LT.3 ) RETURN
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      K = 0
      WRITE (1337,FMT=99005) 'real-space'
      DO L = 1,NSHLR
         WRITE (1337,99006) L,NSR(L),RMR(K+1),(RM(M,K+1),M=1,3)
         DO N = 2,NSR(L)
            WRITE (1337,FMT=99007) (RM(M,K+N),M=1,3)
         END DO
         IF ( L.NE.NSHLR ) WRITE (1337,99008)
         K = K + NSR(L)
      END DO
      WRITE (1337,99009)
      K = 0
      WRITE (1337,FMT=99005) 'reciprocal'
      DO L = 1,NSHLG
         WRITE (1337,99006) L,NSG(L),GNR(K+1),(GN(M,K+1),M=1,3)
         DO N = 2,NSG(L)
            WRITE (1337,FMT=99007) (GN(M,K+N),M=1,3)
         END DO
         IF ( L.NE.NSHLG ) WRITE (1337,99008)
         K = K + NSG(L)
      END DO
      WRITE (1337,99009)
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
99001 FORMAT (10X,'R max =',F9.5,' (a.u.)',/,10X,'G max =',F9.5,
     &        ' (1/a.u.)',/)
99002 FORMAT (10X,'               vectors  shells  max. R ',/,10X,
     &        '               ------------------------------')
99003 FORMAT (10X,A,I7,2X,I6,2X,F9.5)
99004 FORMAT (10X,'               ------------------------------',/)
99005 FORMAT (10X,55('+'),/,18X,'generated ',A,' lattice vectors',/,10X,
     &        55('+'),/,10X,
     &        'shell Nvec    radius          x         y         z',/,
     &        10X,55('-'))
99006 FORMAT (10X,I5,I5,F12.6,2X,3F10.5)
99007 FORMAT (34X,3F10.5)
99008 FORMAT (13X,52('-'))
99009 FORMAT (10X,55('+'),/)
      END
