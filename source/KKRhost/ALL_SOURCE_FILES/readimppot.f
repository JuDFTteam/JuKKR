      SUBROUTINE READIMPPOT(NATOMIMP,INS,IPF,IPFE,IPE,KWS,NSPIN,LPOT,
     +                  IPANIMP,THETASIMP,IRCUTIMP,IRWSIMP,KHFELD,
     +                  HFIELD,VINSIMP,VM2ZIMP,IRMINIMP,RIMP,ZIMP)
c ************************************************************************
c read in impurity potential
c n.h.long, May 2013
c-----------------------------------------------------------------------
C     .. Parameters ..
      IMPLICIT NONE
      include 'inc.p'
      INTEGER NSPIN,NATOMIMP
      INTEGER NPOTD
      PARAMETER (NPOTD=NSPIND*NATOMIMPD)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER LMXSPD
      PARAMETER (LMXSPD= (2*LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT,HFIELD,VBC(2)
      INTEGER INS,IPE,IPF,IPFE,
     +        KHFELD,KWS,LPOT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NATOMIMP),B(NATOMIMP),DRDI(IRMD,NATOMIMP),
     +                 DROR(IRMD,NATOMIMP),ECORE(20,NPOTD),
     +                 RIMP(IRMD,NATOMIMP),RMT(NATOMIMP),
     +                 RMTNEW(NATOMIMP),RWS(NATOMIMP),
     +                 THETASIMP(IRID,NFUND,NATOMIMP),
     +                 VINSIMP(IRMIND:IRMD,LMPOTD,NATOMIMP*NSPIN),
     +                 VM2ZIMP(IRMD,NATOMIMP*NSPIN),ZIMP(NATOMIMP)
      INTEGER IMT(NATOMIMP),IPANIMP(NATOMIMP),
     +        IRCUTIMP(0:IPAND,NATOMIMP),
     +        IRMINIMP(NATOMIMP),IRWSIMP(NATOMIMP),ITITLE(20,NPOTD),
     +        LCORE(20,NPOTD),
     +        NCORE(NPOTD),NFU(NATOMIMP)
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION dummy2(IRMD,NATOMIMP*NSPIN)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,B1,EA,EFNEW,DUMMY
      INTEGER I,IA,ICELL,ICORE,IFUN,IH,IMT1,INEW,IO,IPAN1,IR,IRC1,IRI,
     +        IRMINM,IRMINP,IRNS1P,IRT1P,IRWS1,ISAVE,ISPIN,ISUM,
     +        J,LM,LM1,LMPOT,LMPOTP,
     +        N,NCELL,NFUN,NR
      LOGICAL TEST
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DRN(IRID,NATOMIMP),SCALE(1),U(IRMD),
     +                 XRN(IRID,NATOMIMP)
      INTEGER MESHN(NATOMIMP),NM(IPAND,NATOMIMP),NPAN(NATOMIMP)
C     ..
C     .. External Subroutines ..
      EXTERNAL CALRMT,POTCUT,RINIT,TEST
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ANINT,EXP,LOG,MAX,MOD,REAL,SQRT
C     ..
c-----------------------------------------------------------------------
      WRITE(1337,*) 'in readimppot'
      write(*,*) shape(vinsimp)
      VINSIMP=0d0
c-----------------------------------------------------------------------
c read data from shapefun_imp file
      IF (INS.GT.0) THEN
      OPEN(UNIT=20,FILE='shapefun_imp',FORM='FORMATTED')
      READ (20,*) NCELL
      READ (20,*) SCALE(1)
      DO 30 ICELL = 1,NCELL
       READ (20,FMT=9000) NPAN(ICELL),MESHN(ICELL)
       READ (20,FMT=9000) (NM(IPAN1,ICELL),IPAN1=2,NPAN(ICELL)+1)
       READ (20,FMT=9010) (XRN(IR,ICELL),DRN(IR,ICELL),IR=1,
     +      MESHN(ICELL))
       READ (20,FMT=9000) NFU(ICELL)
       NFUN = NFU(ICELL)

        DO 20 IFUN = 1,NFUN
         READ (20,FMT=9000) LM
         IF (LM.LE.LMXSPD) THEN
          READ (20,FMT=9010) (THETASIMP(N,IFUN,ICELL),N=1,MESHN(ICELL))
         ELSE
          READ (20,FMT=9010) (DUMMY,N=1,MESHN(ICELL))
         END IF
   20     CONTINUE

   30   CONTINUE
      END IF                        ! INS.EQ.1

      DO ICELL=1,NCELL
       IF (INS.NE.0) THEN
        IPANIMP(ICELL) = 1 + NPAN(ICELL) 
       ELSE
        IPANIMP(ICELL)=1
       ENDIF
      ENDDO
c-----------------------------------------------------------------------
c read in impurity potential

      OPEN(UNIT=21,FILE='potential_imp',FORM='FORMATTED')
      LMPOT = (LPOT+1)* (LPOT+1)
      DO 150 IH = 1,NCELL
       DO 140 ISPIN = 1,NSPIN
        I = NSPIN* (IH-1) + ISPIN
        IRCUTIMP(0,IH) = 0

c---> read title of potential card
        READ (21,FMT=9020) (ITITLE(IA,I),IA=1,20)
c
c---  >read muffin-tin radius , lattice constant and new muffin radius
        !READ (21,FMT=9030) RMT(IH),ALAT,RMTNEW(IH)
        READ (21,FMT=*) RMT(IH),ALAT,RMTNEW(IH)
c
c---> read nuclear charge , lmax of the core states ,
c     wigner seitz radius , fermi energy and energy difference
c     between electrostatic zero and muffin tin zero
c
        !READ (21,FMT=9040) ZIMP(IH),RWS(IH),EFNEW,VBC(ISPIN)
        READ (21,FMT=*) ZIMP(IH)
        READ (21,FMT=*) RWS(IH),EFNEW,VBC(ISPIN)
c
c---> read : number of radial mesh points
c     (in case of ws input-potential: last mesh point corresponds
c     to ws-radius, in case of shape-corrected input-potential
c     last mesh point of the exponential mesh corresponds to
c     mt-radius/nevertheless this point is always in the array
c     irws(ih)),number of points for the radial non-muffin-tin
c     mesh  needed for shape functions, the constants a and b
c     for the radial exponential mesh : r(i) = b*(exp(a*(i-1))-1)
c     the no. of different core states and some other stuff
c
        READ (21,FMT=9050) IRWSIMP(IH)
        !READ (21,FMT=9051) A(IH),B(IH),NCORE(I),INEW
        READ (21,FMT=*) A(IH),B(IH)
        READ (21,FMT=*) NCORE(I),INEW
        NR = IRWSIMP(IH)
c---> read the different core states : l and energy
c
        IF (NCORE(I).GE.1) READ (21,FMT=9070) (LCORE(ICORE,I),
     +      ECORE(ICORE,I),ICORE=1,NCORE(I))
c
         IF (INS.LT.1) THEN
c
c--->  read radial mesh points, its derivative, the spherically averaged
c      charge density and the input potential without the nuclear pot.
c
          IF (INEW.EQ.0) THEN
           READ (21,FMT=9060) (RIMP(IR,IH),DRDI(IR,IH),VM2ZIMP(IR,I),
     +           IR=1,NR)
          ELSE
           READ (21,FMT=*) (VM2ZIMP(IR,I),IR=1,NR)
          END IF

         ELSE                    ! (INS.LT.1)
c
c--->  read full potential - the non spherical contribution from irmin
c      to irt - remember that the lm = 1 contribution is multiplied by
c      1/sqrt(4 pi)
c
          READ (21,FMT=9090) IRT1P,IRNS1P,LMPOTP,ISAVE
          IRMINP = IRT1P - IRNS1P
          IRMINM = MAX(IRMINP,IRMIND)
          READ (21,FMT=9100) (VM2ZIMP(IR,I),IR=1,NR)
          IF (LMPOTP.GT.1) THEN
           LM1 = 2
           DO 50 LM = 2,LMPOTP
            IF (LM1.NE.1) THEN
             IF (ISAVE.EQ.1) THEN
              READ (21,FMT=9090) LM1
             ELSE
              LM1 = LM
             END IF
             IF (LM1.GT.1) THEN
              READ (21,FMT=9100) (U(IR),IR=IRMINP,NR)
              IF (LM1.LE.LMPOT) THEN
               DO 40 IR = IRMINM,NR
                VINSIMP(IR,LM1,I) = U(IR)
   40          CONTINUE
              END IF
             END IF
            END IF
   50      CONTINUE
          END IF
         END IF                  ! (INS.LT.1)
         IRWS1 = IRWSIMP(IH)

c---> redefine new mt-radius in case of shape corrections
c
         IF (INS.NE.0) THEN
          RMTNEW(IH) = SCALE(1)*ALAT*XRN(1,IH)
          IMT1 = ANINT(LOG(RMTNEW(IH)/B(IH)+1.0D0)/A(IH)) + 1
c
c---> for proper core treatment imt must be odd
c     shift potential by one mesh point if imt is even
c
          IF (MOD(IMT1,2).EQ.0) THEN
           IMT1 = IMT1 + 1
           DO 60 IR = IMT1,2,-1
            VM2ZIMP(IR,I) = VM2ZIMP(IR-1,I)
   60      CONTINUE
          END IF

          IMT(IH) = IMT1
          B(IH) = RMTNEW(IH)/ (EXP(A(IH)*REAL(IMT1-1))-1.0D0)
         END IF               ! (INS.NE.0)

c---> generate radial mesh - potential only is stored in potential card
c     INEW = 1
c     p. zahn, jan. 99
c
         A1 = A(IH)
         B1 = B(IH)
         RIMP(1,IH) = 0.0D0
         DRDI(1,IH) = A1*B1
         DO 70 IR = 2,IRWS1
          EA = EXP(A1*REAL(IR-1))
          RIMP(IR,IH) = B1* (EA-1.0D0)
          DRDI(IR,IH) = A1*B1*EA
          DROR(IR,IH) = A1/ (1.0D0-1.0D0/EA)
   70    CONTINUE
c
c---> fill cell-type depending mesh points in the non-muffin-tin-region
c
          IF (INS.NE.0) THEN
           DO 80 IRI = 1,MESHN(IH)
            IR = IRI + IMT1
            RIMP(IR,IH) = SCALE(1)*ALAT*XRN(IRI,IH)
            DRDI(IR,IH) = SCALE(1)*ALAT*DRN(IRI,IH)
            DROR(IR,IH) = DRDI(IR,IH)/RIMP(IR,IH)
   80      CONTINUE
          END IF

          RWS(IH) = RIMP(IRWS1,IH)
c
c---> kshape.eq.0 : calculate new rmt adapted to exp. mesh
c
          CALL CALRMT(IPF,IPFE,IPE,IMT(IH),ZIMP(IH),RMT(IH),RWS(IH),
     +                RMTNEW(IH),ALAT,DRDI(1,IH),A(IH),B(IH),IRWS1,
     +                RIMP(1,IH),IO,INS)

          IF (INS.GT.0) THEN
           IRCUTIMP(1,IH) = IMT(IH)
           ISUM = IMT(IH)           
            DO 90 IPAN1 = 2,IPANIMP(IH)
             ISUM = ISUM + NM(IPAN1,IH)
             IRCUTIMP(IPAN1,IH) = ISUM              
   90       CONTINUE
            NR = ISUM
          ELSE                    ! INS.EQ.0
            NR = IRWSIMP(IH)
            IF (KWS.GE.1) THEN
             IRCUTIMP(1,IH) = IRWS1
            ELSE
             IRCUTIMP(1,IH) = IMT(IH)
            END IF
           END IF                  ! INS.GT.0
c
            
c
c---> fill array irmin in case of full potential
c
          
          IF (INS.NE.0) IRMINIMP(IH) = NR - IRNS1P

c---> cut input potential at rmt if given only at exponential mesh
c
          IF (INS.GT.1) THEN
           IMT1 = IMT(IH)
           IRC1 = IRCUTIMP(IPANIMP(IH),IH)
           CALL POTCUT(IMT1,IRC1,INS,LMPOT,RIMP(1,IH),VM2ZIMP(1,I),
     +                 dummy2,VINSIMP(IRMIND,1,I),ZIMP(IH),IRMD,IRMIND)
          END IF
c

c
          IF (INS.EQ.0 .AND. KWS.EQ.0) THEN
c
c---> in case of a mt calculation cut potential at mt radius
c
           IMT1 = IMT(IH)
           IRWS1 = IRWSIMP(IH)
           CALL POTCUT(IMT1,IRWS1,INS,LMPOT,RIMP(1,IH),VM2ZIMP(1,I),
     +                 dummy2,VINSIMP(IRMIND,1,I),ZIMP(IH),IRMD,IRMIND)

          END IF                    ! INS.EQ.0 .AND. KWS.EQ.0
c
c--->       maybe apply a magnetic field
          IF (KHFELD.EQ.1 ) THEN
           WRITE(1337,*) 'ATOM',IH,'SPIN',ISPIN,'SHIFTED BY',
     +           -REAL(2*ISPIN-3)*HFIELD
           DO 130 J = 1,IRCUTIMP(IPANIMP(IH),IH)
            VM2ZIMP(J,I) = VM2ZIMP(J,I) - REAL(2*ISPIN-3)*HFIELD
  130      CONTINUE
          END IF

  140    CONTINUE                    ! ISPIN = 1,NSPIN
  150  CONTINUE                      ! IH = 1,NCELL
      CLOSE(20)
      CLOSE(21)

      RETURN


 9000 FORMAT (16i5)
 9010 FORMAT (4d20.12)
 9020 FORMAT (20a4)
 9030 FORMAT (3f12.8)
 9040 FORMAT (f10.5,/,f10.5,2f15.10)
 9050 FORMAT (i4)
 9051 FORMAT (2d15.8,/,2i2)
 9060 FORMAT (1p,2d15.6,1p,d15.8)
 9070 FORMAT (i5,1p,d20.11)
 9090 FORMAT (10i5)
 9100 FORMAT (1p,4d20.13)
      END                      
