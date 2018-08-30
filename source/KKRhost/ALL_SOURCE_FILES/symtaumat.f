C*==symtaumat.f    processed by SPAG 6.05Rc at 15:50 on 10 Dec 2002
      SUBROUTINE SYMTAUMAT(ROTNAME,ROTMAT,DROT,NSYM,ISYMINDEX,
     &                     SYMUNITARY,NQMAX,NKMMAX,NQ,NL,KREL,IPRINT,
     &                     NSYMAXD)
C   ********************************************************************
C   *                                                                  *
C   *  Find the symmetry matrices DROT that act on t, tau, ....        *
C   *  KREL=0: for real spherical harmonics                            *
C   *  KREL=1: for relativistic represntation                          *
C   *                                                                  *
C   *  The NSYM allowed symmetry operations are indicated by ISYMINDEX *
C   *  in the table  ROTMAT. For KREL=1, SYMUNITARY=T/F indicates a    *
C   *  unitary/antiunitary symmetry operation.                         *
C   *                                                                  *
C   *  The routine determines first the Euler angles correponding      *
C   *  to a symmetry operation. Reflections are decomposed into        *
C   *  inversion + rotation for this reason.                           *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C PARAMETER definitions
C
      DOUBLE COMPLEX CI,C1,C0
      PARAMETER (CI=(0.0D0,1.0D0),C1=(1.0D0,0.0D0),C0=(0.0D0,0.0D0))
C
C Dummy arguments
C
      INTEGER IPRINT,KREL,NKMMAX,NL,NQ,NQMAX,NSYM,NSYMAXD
      DOUBLE COMPLEX DROT(NKMMAX,NKMMAX,48)
      INTEGER ISYMINDEX(NSYMAXD)
      DOUBLE PRECISION ROTMAT(64,3,3)
      CHARACTER (len=10) :: ROTNAME(64)
      LOGICAL SYMUNITARY(48)
C
C Local variables
C
      double precision  A,B,CO1,CO2,CO3,DET,FACT(0:100),PI,RJ,RMJ,
     &        SI1,SI2,SI3,SK,SYMEULANG(3,48),TET1,TET2,TET3
      LOGICAL CHECKRMAT
      DOUBLE PRECISION DBLE
      double precision  DDET33
      DOUBLE COMPLEX DINV(NKMMAX,NKMMAX),DTIM(NKMMAX,NKMMAX),
     &           RC(NKMMAX,NKMMAX),W1(NKMMAX,NKMMAX),W2(NKMMAX,NKMMAX)
      LOGICAL EQUAL
      INTEGER I,I1,I2,IND0Q(NQMAX),INVFLAG(48),IQ,IREL,IRELEFF,ISYM,
     &        ITOP,J,K,L,LOOP,M,N,NK,NKEFF,NKM,NLM,NOK,NS
      INTEGER NINT
      DOUBLE PRECISION RMAT(3,3)
      double precision  W
C
      EQUAL(A,B) = (DABS(A-B).LT.1D-7)
C
      WRITE (1337,99001)
C
      PI = 4D0*ATAN(1D0)
C
      IREL = KREL*3
      NK = (1-KREL)*NL + KREL*(2*NL-1)
      NKM = (1+KREL)*NL**2
C
C-----------------------------------------------------------------------
      FACT(0) = 1.0D0
      DO I = 1,100
         FACT(I) = FACT(I-1)*DBLE(I)
      END DO
C-----------------------------------------------------------------------
C
      IND0Q(1) = 0
      DO IQ = 2,NQ
         IND0Q(IQ) = IND0Q(IQ-1) + NKM
      END DO
C
C ----------------------------------------------------------------------
C    RC  transforms from  REAL to  COMPLEX (L,M,S) - representation
C                 |LC> = sum[LR] |LR> * RC(LR,LC)
C ----------------------------------------------------------------------
      IF ( KREL.EQ.0 ) THEN
         NLM = NKM
C
         CALL CINIT(NKMMAX*NKMMAX,RC)
C
         W = 1.0D0/SQRT(2.0D0)
C
         DO L = 0,(NL-1)
            DO M = -L,L
               I = L*(L+1) + M + 1
               J = L*(L+1) - M + 1
C
               IF ( M.LT.0 ) THEN
                  RC(I,I) = -CI*W
                  RC(J,I) = W
               END IF
               IF ( M.EQ.0 ) THEN
                  RC(I,I) = C1
               END IF
               IF ( M.GT.0 ) THEN
                  RC(I,I) = W*(-1.0D0)**M
                  RC(J,I) = CI*W*(-1.0D0)**M
               END IF
            END DO
         END DO
      END IF
C
C=======================================================================
C     The routine determines first the Euler angles correponding
C     to a symmetry operation. Reflections are decomposed into
C     inversion + rotation for this reason.
C=======================================================================
C
      DO ISYM = 1,NSYM
C
         DO I1 = 1,3
            DO I2 = 1,3
               RMAT(I1,I2) = ROTMAT(ISYMINDEX(ISYM),I1,I2)
            END DO
         END DO
C
         DET = DDET33(RMAT)
C
         INVFLAG(ISYM) = 0
         IF ( DET.LT.0D0 ) THEN
            CALL DSCAL(9,-1.0D0,RMAT,1)
            INVFLAG(ISYM) = 1
         END IF
C
C----------------------------------------------------------------------
         CO2 = RMAT(3,3)
         TET2 = ACOS(CO2)
         LOOP = 0
 50      CONTINUE
         IF ( LOOP.EQ.1 ) TET2 = -TET2
         SI2 = SIN(TET2)
C
         IF ( EQUAL(CO2,1.0D0) ) THEN
            TET1 = ACOS(RMAT(1,1))
            IF ( .NOT.EQUAL(RMAT(1,2),SIN(TET1)) ) THEN
               TET1 = -TET1
               IF ( .NOT.EQUAL(RMAT(1,2),SIN(TET1)) ) WRITE (1337,*)
     &               '>>>>>>>>>>>>>>> STRANGE 1'
            END IF
            TET2 = 0D0
            TET3 = 0D0
         ELSE IF ( EQUAL(CO2,-1D0) ) THEN
            TET1 = ACOS(-RMAT(1,1))
            IF ( .NOT.EQUAL(RMAT(1,2),-SIN(TET1)) ) THEN
               TET1 = -TET1
               IF ( .NOT.EQUAL(RMAT(1,2),-SIN(TET1)) ) WRITE (1337,*)
     &               '>>>>>>>>>>>>>>> STRANGE 2'
            END IF
            TET2 = PI
            TET3 = 0D0
         ELSE
            TET1 = ACOS(RMAT(3,1)/SI2)
            IF ( .NOT.EQUAL(RMAT(3,2),SI2*SIN(TET1)) ) THEN
               TET1 = -TET1
               IF ( .NOT.EQUAL(RMAT(3,2),SI2*SIN(TET1)) ) WRITE (1337,*)
     &               '>>>>>>>>>>>>>>> STRANGE 3'
            END IF
C
            TET3 = ACOS(-RMAT(1,3)/SI2)
            IF ( .NOT.EQUAL(RMAT(2,3),SI2*SIN(TET3)) ) THEN
               TET3 = -TET3
               IF ( .NOT.EQUAL(RMAT(2,3),SI2*SIN(TET3)) ) WRITE (1337,*)
     &               '>>>>>>>>>>>>>>> STRANGE 4'
            END IF
C
         END IF
C
         CO1 = COS(TET1)
         SI1 = SIN(TET1)
         CO2 = COS(TET2)
         SI2 = SIN(TET2)
         CO3 = COS(TET3)
         SI3 = SIN(TET3)
C
         NOK = 0
         DO I1 = 1,3
            DO I2 = 1,3
               IF ( CHECKRMAT(RMAT,CO1,SI1,CO2,SI2,CO3,SI3,I1,I2) ) THEN
                  NOK = NOK + 1
               ELSE IF ( LOOP.LT.1 ) THEN
                  LOOP = LOOP + 1
                  GOTO 50
               END IF
            END DO
         END DO
C
         SYMEULANG(1,ISYM) = TET1*(180D0/PI)
         SYMEULANG(2,ISYM) = TET2*(180D0/PI)
         SYMEULANG(3,ISYM) = TET3*(180D0/PI)
C
         IF ( NOK.NE.9 ) WRITE (1337,99009) NOK
         WRITE (1337,99008) ISYM,ROTNAME(ISYMINDEX(ISYM)),INVFLAG(ISYM),
     &                   (SYMEULANG(I,ISYM),I=1,3),SYMUNITARY(ISYM)
C
      END DO
      WRITE(1337,'(8X,57(1H-),/)')
C
C-----------------------------------------------------------------------
C                    initialize all rotation matrices
C-----------------------------------------------------------------------
C
      CALL CINIT(NKMMAX*NKMMAX*NSYM,DROT)
C
C-----------------------------------------------------------------------
C                       create rotation matrices
C-----------------------------------------------------------------------
C
      IF ( IREL.LE.2 ) THEN
         IRELEFF = 0
         NKEFF = NL
      ELSE
         IRELEFF = 3
         NKEFF = NK
      END IF
C
      DO ISYM = 1,NSYM
C
         CALL CALCROTMAT(NKEFF,IRELEFF,SYMEULANG(1,ISYM),
     &                   SYMEULANG(2,ISYM),SYMEULANG(3,ISYM),
     &                   DROT(1,1,ISYM),FACT,NKMMAX)
C
      END DO
C-----------------------------------------------------------------------
C                     create matrix for inversion
C-----------------------------------------------------------------------
      CALL CINIT(NKMMAX*NKMMAX,DINV)
C
      I = 0
      IF ( IREL.GT.2 ) THEN
         NS = 2
      ELSE
         NS = 1
      END IF
      DO L = 0,(NL-1)
         DO M = 1,NS*(2*L+1)
            I = I + 1
            DINV(I,I) = (-1.0D0)**L
         END DO
      END DO
      ITOP = I
C
C-----------------------------------------------------------------------
C                         include inversion
C-----------------------------------------------------------------------
      DO ISYM = 1,NSYM
         IF ( INVFLAG(ISYM).NE.0 ) THEN
C
            CALL ZGEMM('N','N',NKM,NKM,NKM,C1,DROT(1,1,ISYM),NKMMAX,
     &                 DINV,NKMMAX,C0,W2,NKMMAX)
C
            DO J = 1,NKM
               CALL ZCOPY(NKM,W2(1,J),1,DROT(1,J,ISYM),1)
            END DO
         END IF
      END DO
C
C-----------------------------------------------------------------------
C            add second spin-diagonal block for  IREL=2
C            spin off-diagonal blocks have been initialized before
C-----------------------------------------------------------------------
      IF ( IREL.EQ.2 ) THEN
         NLM = NKM/2
         IF ( ITOP.NE.NLM ) CALL ERRORTRAP('SYMTAUMAT',11,1)
         DO ISYM = 1,NSYM
C
            DO J = 1,NLM
               CALL ZCOPY(NLM,DROT(1,J,ISYM),1,DROT(NLM+1,NLM+J,ISYM),1)
            END DO
         END DO
      END IF
C-----------------------------------------------------------------------
C            transform to real spherical representation for  KREL=0
C-----------------------------------------------------------------------
      N = NKM
      M = NKMMAX
      IF ( KREL.EQ.0 ) THEN
         DO ISYM = 1,NSYM
            CALL ZGEMM('N','N',N,N,N,C1,RC,M,DROT(1,1,ISYM),M,C0,W1,M)
            CALL ZGEMM('N','C',N,N,N,C1,W1,M,RC,M,C0,DROT(1,1,ISYM),M)
         END DO
      END IF
C-----------------------------------------------------------------------
C                     create matrix for time reversal
C-----------------------------------------------------------------------
      IF ( IREL.GT.1 ) THEN
C
         CALL CINIT(NKMMAX*NKMMAX,DTIM)
C
         I = 0
         DO K = 1,NK
            L = K/2
            IF ( L*2.EQ.K ) THEN
               SK = -1D0
            ELSE
               SK = +1D0
            END IF
            RJ = L + SK*0.5D0
            DO RMJ = -RJ, + RJ
               I1 = NINT(2*L*(RJ+0.5D0)+RJ+RMJ+1)
               I2 = NINT(2*L*(RJ+0.5D0)+RJ-RMJ+1)
               DTIM(I1,I2) = SK*(-1)**NINT(RMJ+0.5D0)
            END DO
         END DO
         IF ( IPRINT.GT.0 ) THEN
            CALL CMATSTR('Inversion     MATRIX',20,DINV,NKM,NKMMAX,3,3,
     &                   0,1D-8,6)
            CALL CMATSTR('Time reversal MATRIX',20,DTIM,NKM,NKMMAX,3,3,
     &                   0,1D-8,6)
         END IF
C
      END IF
C=======================================================================
C            set up of transformation matrices completed
C=======================================================================
C
C=======================================================================
C   include time reversal operation for anti-unitary transformations
C=======================================================================
      DO ISYM = 1,NSYM
         IF ( .NOT.SYMUNITARY(ISYM) ) THEN
            IF ( IREL.EQ.2 ) CALL ERRORTRAP('SYMTAUMAT',14,1)
C
            CALL ZGEMM('N','N',NKM,NKM,NKM,C1,DROT(1,1,ISYM),NKMMAX,
     &                 DTIM,NKMMAX,C0,W2,NKMMAX)
            DO J = 1,NKM
               CALL ZCOPY(NKM,W2(1,J),1,DROT(1,J,ISYM),1)
            END DO
         END IF
      END DO
C
C-----------------------------------------------------------------------
C for testing
C
Cccc      write (6,*) ' NUMBER OF SYMMETRIES : ', NSYM
Cccc      
Cccc      do isym = 1,nsym
Cccc         write(6,*) ' ISYM = ',isym
Cccc         call cmatstr('DROT',4,drot(1,1,isym),nkm,nkmmax,krel*3,krel*3,
Cccc     &        0,1d-12,6)
Cccc         write(6,*)
Cccc      end do
C     
C-----------------------------------------------------------------------
C
      IF ( IPRINT.EQ.0 ) RETURN
C
C=======================================================================
C       find the structure of the site-diagonal TAU - matrices  TAUQ
C=======================================================================
C
      CALL TAUSTRUCT(DROT,NSYM,SYMUNITARY,NKM,NQ,NQMAX,NKMMAX,IPRINT,
     &               IREL)
C
      RETURN
99001 FORMAT (5X,'<SYMTAUMAT> : rotation matrices acting on t/G/tau',//,
     &        8X,57(1H-),/,8X,
     &     'ISYM            INV          Euler angles      Unitarity',/,
     &        8X,57(1H-))
99008 FORMAT (8X,I2,3X,A,I3,3F10.5,3X,L1)
99009 FORMAT (50('>'),' trouble in <SYMTAUMAT>',I3,F10.5)
      END
