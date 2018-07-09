      SUBROUTINE HFFCORE(RNUC,JTOP,KAP1,KAP2,NSOL,MJ,GC,FC,NRC,SHF,S,
     &                   NMEMAX,NKMMAX,R,DRDI,SDIA,SMDIA,SOFF,SMOFF,
     &                   QDIA,QOFF,QMDIA,QMOFF,NUCLEUS,JLIM)
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Calculates matrix elements of several hyperfine interaction
C     connected quantities in the core.
C     All the related arrays have a counting index as
C     the last index of the array indicates the corresponding physical
C     property.
C     Index-list
C     1      electron-Spin-electron-Spin Hyperfine field
C     2      nuclear-spin-electron-orbit hyperfine field
C     3      electron-spin-nulceus-spin-contact hyperfine field
C     4      expectation value of (1/r)^3
C     5      Total Hyperfine Field (see Rose (1961))
C     called by core
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
C
C
C PARAMETER definitions
C
      DOUBLE PRECISION MB,A0,F1,F2
C
C     BOHR-MAGNETON       IN ERG/GAUSS
C     CONVERSION FACTOR FOR HYPERFINE FIELDS FROM A.U. TO GAUSS
C     ELECTRON CHARGE     IN ESU
C
C
      PARAMETER (MB=9.274078D-21,A0=0.52917706D-08,F1=1.0D0,
     &           F2=2.0D0*MB/(A0*A0*A0))
C
C Dummy arguments
C
      INTEGER JLIM,JTOP,KAP1,KAP2,NKMMAX,NMEMAX,NRC,NSOL,NUCLEUS,S
      DOUBLE PRECISION MJ,RNUC
      DOUBLE PRECISION DRDI(NRC),FC(2,2,NRC),GC(2,2,NRC),QDIA(NKMMAX),
     &       QMDIA(NKMMAX),QMOFF(NKMMAX),QOFF(NKMMAX),R(NRC),
     &       SDIA(NKMMAX),SHF(2,2,NMEMAX),SMDIA(NKMMAX),SMOFF(NKMMAX),
     &       SOFF(NKMMAX)
C
C Local variables
C
      DOUBLE PRECISION AME(2,2),CFF(2,2),CFG(2,2),CGF(2,2),CGG(2,2),
     &       CQF(2,2),
     &       CQG(2,2),CSF(2,2),CSG(2,2),DOVR(NRC),DROVRN(NRC),
     &       DROVRN1(NRC),F(NRC,2),FF(2,2),FF1(2,2),FF2(2,2),FG(2,2),
     &       FG1(2,2),FG2(2,2),G(NRC,2),GF(2,2),GF1(2,2),GF2(2,2),
     &       GG(2,2),GG1(2,2),GG2(2,2)
      DOUBLE PRECISION DBLE,DSQRT
      INTEGER I,IKM1,IKM2,J,K,K1,K2,N
      INTEGER IKAPMUE
      INTEGER NINT
C
      IF ( KAP2.EQ.0 ) KAP2 = KAP1
C
      CALL RINIT(4,GG)
      CALL RINIT(4,FF)
      CALL RINIT(4,GG1)
      CALL RINIT(4,FF1)
      CALL RINIT(4,GG2)
      CALL RINIT(4,FF2)
      CALL RINIT(4,GF)
      CALL RINIT(4,FG)
      CALL RINIT(4,GF1)
      CALL RINIT(4,FG1)
      CALL RINIT(4,GF2)
      CALL RINIT(4,FG2)
      CALL RINIT(2*NRC,G)
      CALL RINIT(2*NRC,F)
C
      DO K = 1,2
         DO N = 1,JTOP
            G(N,K) = GC(K,S,N)
            F(N,K) = FC(K,S,N)
         END DO
      END DO
C     prepare meshes for finite nucleus calculation
      DO I = 1,JTOP
         DOVR(I) = DRDI(I)/R(I)
         IF ( NUCLEUS.NE.0 ) THEN
            DROVRN1(I) = (R(I)/RNUC)**3*DRDI(I)
            DROVRN(I) = DROVRN1(I)/R(I)
         END IF
      END DO
      IKM1 = IKAPMUE(KAP1,NINT(MJ-0.5D0))
      IKM2 = IKAPMUE(KAP2,NINT(MJ-0.5D0))
C     angular hyperfine matrix elements   see e.g.  E.M.Rose
C     the factor  i  has been omitted
      CALL RINIT(4,AME)
      AME(1,1) = 4.0D0*KAP1*MJ/(4.0D0*KAP1*KAP1-1.0D0)
      IF ( NSOL.EQ.2 ) THEN
         AME(2,2) = 4.0D0*KAP2*MJ/(4.0D0*KAP2*KAP2-1.0D0)
         AME(1,2) = DSQRT(0.25D0-(MJ/DBLE(KAP1-KAP2))**2)
         AME(2,1) = AME(1,2)
      END IF
C     coefficients for the spin-dipolar matrix elements
      CALL RINIT(4,CSF)
      CALL RINIT(4,CSG)
      CSG(1,1) = SDIA(IKM1)
      CSF(1,1) = SMDIA(IKM1)
      IF ( NSOL.EQ.2 ) THEN
         CSG(2,2) = SDIA(IKM2)
         CSG(1,2) = SOFF(IKM1)
         CSG(2,1) = CSG(1,2)
         CSF(2,2) = SMDIA(IKM2)
         CSF(1,2) = SMOFF(IKM1)
         CSF(2,1) = SMOFF(IKM1)
      END IF
C     COEFFICIENTS FOR THE QUADRUPOLAR MATRIX ELEMENTS
      CQG(1,1) = QDIA(IKM1)
      CQG(2,2) = QDIA(IKM2)
      CQG(1,2) = QOFF(IKM1)
      CQG(2,1) = CQG(1,2)
      CALL RINIT(4,CQF)
      CQF(1,1) = QMDIA(IKM1)
      CQF(2,2) = QMDIA(IKM2)
      CQF(1,2) = QMOFF(IKM1)
      CQF(2,1) = CQF(1,2)
C     coefficients to calculate the spin-spin field
      CALL RINIT(4,CGG)
      CALL RINIT(4,CGF)
      CGG(1,1) = -MJ/(KAP1+0.5D0)
      CGF(1,1) = -MJ/(-KAP1+0.5D0)
      IF ( NSOL.EQ.2 ) THEN
         CGG(1,2) = -DSQRT(1.0D0-(MJ/(KAP1+0.5D0))**2)
         CGG(2,1) = CGG(1,2)
         CGG(2,2) = -MJ/(KAP2+0.5D0)
         CGF(2,2) = -MJ/(-KAP2+0.5D0)
C     CGF(1,2) = -DSQRT( 1.0D0 - (MJ/(- KAP1+0.5D0))**2 )
C     CGF(2,1) = CGF(1,2)
      END IF
C     coefficients to calculate the orbital field
      CALL RINIT(4,CFG)
      CALL RINIT(4,CFF)
      CFG(1,1) = MJ*(KAP1+1.0D0)/(KAP1+0.5D0)
      CFF(1,1) = MJ*(-KAP1+1.0D0)/(-KAP1+0.5D0)
      IF ( NSOL.EQ.2 ) THEN
         CFG(2,2) = MJ*(KAP2+1.0D0)/(KAP2+0.5D0)
         CFG(1,2) = 0.5D0*DSQRT(1.0D0-(MJ/(KAP1+0.5D0))**2)
         CFG(2,1) = CFG(1,2)
         CFF(2,2) = MJ*(-KAP2+1.0D0)/(-KAP2+0.5D0)
C     CFF(1,2) = 0.5D0 * DSQRT( 1.0D0 - (MJ/(- KAP1 + 0.5D0))**2 )
C     CFF(2,1) = CFF(1,2)
      END IF
C Calculates integrals from 0.0 to jtop
      CALL HFFINT(GG,G,G,DOVR,R,0.0D0,NSOL,JTOP,NRC)
      CALL HFFINT(FF,F,F,DOVR,R,0.0D0,NSOL,JTOP,NRC)
      CALL HFFINT(GF,G,F,DRDI,R,0.0D0,NSOL,JTOP,NRC)
      CALL HFFINT(FG,F,G,DRDI,R,0.0D0,NSOL,JTOP,NRC)
      CALL RSUMUPINT(SHF(1,1,5),F1,GG,CQG,F1,FF,CQF,NSOL)
      IF ( NUCLEUS.NE.0 ) THEN
C     calculates integrals inside nucleus at RNUC in order to get
C     contribution outside the nucleus
         CALL HFFINT(GG1,G,G,DOVR,R,RNUC,NSOL,JLIM,NRC)
         CALL HFFINT(FF1,F,F,DOVR,R,RNUC,NSOL,JLIM,NRC)
         CALL HFFINT(GF1,G,F,DRDI,R,RNUC,NSOL,JLIM,NRC)
         CALL HFFINT(FG1,F,G,DRDI,R,RNUC,NSOL,JLIM,NRC)
C     calculates contribution from RNUC to jtop
         DO I = 1,NSOL
            DO J = 1,NSOL
               GG(I,J) = GG(I,J) - GG1(I,J)
               FF(I,J) = FF(I,J) - FF1(I,J)
               GF(I,J) = GF(I,J) - GF1(I,J)
               FG(I,J) = FG(I,J) - FG1(I,J)
            END DO
         END DO
      END IF                    !end of nucleus.eq.0
C     calculates B_sp which is zero inside the nucleus
      CALL RSUMUPINT(SHF(1,1,1),F1,GG,CSG,-F1,FF,CSF,NSOL)
C     calculates hyperfine integrals from 0.0 to RNUC which are added
C     external integrals
      IF ( NUCLEUS.NE.0 ) THEN
         CALL HFFINT(GG2,G,G,DROVRN,R,RNUC,NSOL,JLIM,NRC)
         CALL HFFINT(FF2,F,F,DROVRN,R,RNUC,NSOL,JLIM,NRC)
         CALL HFFINT(GF2,G,F,DROVRN1,R,RNUC,NSOL,JLIM,NRC)
         CALL HFFINT(FG2,F,G,DROVRN1,R,RNUC,NSOL,JLIM,NRC)
         DO I = 1,NSOL
            DO J = 1,NSOL
               GG(I,J) = GG(I,J) + GG2(I,J)
               FF(I,J) = FF(I,J) + FF2(I,J)
               GF(I,J) = GF(I,J) + GF2(I,J)
               FG(I,J) = FG(I,J) + FG2(I,J)
            END DO
         END DO
      END IF
C     calculates B_nseo and B_tot
      CALL RSUMUPINT(SHF(1,1,2),F2,GG,CFG,-F2,FF,CFF,NSOL)
C      CALL RSUMUPINT(SHF(1,1,5),CAUTOG,GF,AME,CAUTOG,FG,AME,NSOL)
C     modifications for B_ssc which is zero outside the nucleus
      IF ( NUCLEUS.NE.0 ) THEN
         DO I = 1,NSOL
            DO J = 1,NSOL
               GG(I,J) = GG2(I,J)
               FF(I,J) = FF2(I,J)
            END DO
         END DO
      END IF
C     calculates B_ssc
      CALL RSUMUPINT(SHF(1,1,3),F2,GG,CGG,-F2,FF,CGF,NSOL)
C     for testing purposes write in 4 the sum of 1,2,3
      DO K1 = 1,NSOL
         DO K2 = 1,NSOL
            SHF(K1,K2,4) = SHF(K1,K2,1) + SHF(K1,K2,2) + SHF(K1,K2,3)
         END DO
      END DO
C
      END
C
      SUBROUTINE RSUMUPINT(SUM,VG,G,WG,VF,F,WF,N)
      IMPLICIT NONE
C
C
C Dummy arguments
C
      INTEGER N
      DOUBLE PRECISION VF,VG
      DOUBLE PRECISION F(N,N),G(N,N),SUM(N,N),WF(N,N),WG(N,N)
C
C Local variables
C
      INTEGER I,J
C
C
      DO J = 1,N
         DO I = 1,N
            SUM(I,J) = VG*G(I,J)*WG(I,J) + VF*F(I,J)*WF(I,J)
         END DO
      END DO
      END
C
      SUBROUTINE HFFINT(GG,GA,GB,DR,R,RNUC,NSOL,JTOP,NRC)
C     Calculates Hyperfine integrals, extrapolates to zero and
C     intrapolates to exact nuclear radius RNUC
      IMPLICIT NONE
C
C
C Dummy arguments
C
      INTEGER JTOP,NRC,NSOL
      DOUBLE PRECISION RNUC
      DOUBLE PRECISION DR(NRC),GA(NRC,2),GB(NRC,2),GG(2,2),R(NRC)
C
C Local variables
C
      INTEGER I,K1,K2
      DOUBLE PRECISION X(5),Y(5),YI(NRC),ZI(NRC)
      DOUBLE PRECISION YLAG
C
C
      DO K1 = 1,NSOL
         DO K2 = 1,NSOL
            DO I = 1,JTOP
               YI(I) = GA(I,K1)*GB(I,K2)*DR(I)
            END DO
            CALL RINT4PTS(YI,JTOP,ZI)
C     Intrapolation
            IF ( RNUC.NE.0.0D0 ) THEN
               DO I = 1,5
                  X(I) = R(JTOP-5+I)
                  Y(I) = ZI(JTOP-5+I)
               END DO
               ZI(JTOP) = YLAG(RNUC,X,Y,0,4,5)
            END IF
C     Extrapolation
            X(1) = 1.0D0
            X(2) = 6.0D0
            X(3) = 11.0D0
            Y(1) = ZI(JTOP) - ZI(1)
            Y(2) = ZI(JTOP) - ZI(5)
            Y(3) = ZI(JTOP) - ZI(9)
            GG(K1,K2) = YLAG(0.0D0,X,Y,0,2,3)
         END DO
      END DO
      END
C
