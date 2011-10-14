      SUBROUTINE FORCXC(FLM,FLMC,LPOT,NSPIN,IATYP,RHOC,V,R,ALAT,
     +                  DRDI,IRWS,ZAT,
     >                  LMPIC,MYLRANK,
     >                  LGROUP,LCOMM,LSIZE,
C                       new input parameters after inc.p replace
     &                  naez, irmd, prod_lmpid_smpid_empid)
C
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     calculates the force on nucleus m
c     from a given non spherical charge density at the nucleus site r
c     with core correction(exchange contribution)
 
c-----------------------------------------------------------------------

      include 'mpif.h'

      INTEGER naez
      INTEGER irmd
C     LMPID * SMPID * EMPID
      INTEGER prod_lmpid_smpid_empid
C
C      INTEGER LMPOTD
C      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT
      INTEGER LPOT,NSPIN
C     ..
C     .. Array Arguments ..
C      DOUBLE PRECISION DRDI(IRMD,*),FLM(-1:1,*),FLMC(-1:1,*),R(IRMD,*),
C     &       RHOC(IRMD,*),ZAT(NAEZD),
C     &       V(IRMD,LMPOTD,2)

      DOUBLE PRECISION DRDI(IRMD,*),FLM(-1:1,*),FLMC(-1:1,*),R(IRMD,*),
     &       RHOC(IRMD,*),ZAT(NAEZ),
     &       V(IRMD,(LPOT+1)**2,2)

      INTEGER IRWS(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DV,DVOL,FAC,PI,RWS,TRP,VINT1,VOL
      INTEGER I,IATYP,IPER,IPOT,IREP,IRWS1,ISPIN,J,LM,M
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION F(3,NAEZ),
     +                 FALL(3,NAEZ),
     +                 FLMH(-1:1,NAEZ),
     +                 FLMXC(-1:1,NAEZ),V1(IRMD)
C     ..
C     .. External Subroutines ..
      EXTERNAL SIMP3
C     ..
C     .. Save statement ..
      SAVE PI
C     ..
C     .. MPI variables ..
      INTEGER IERR

C     .. L-MPI ..
      INTEGER      MYLRANK(prod_lmpid_smpid_empid),
     +             LCOMM(prod_lmpid_smpid_empid),
     +             LGROUP(prod_lmpid_smpid_empid),
     +             LSIZE(prod_lmpid_smpid_empid),
     +             LMPI,LMPIC
C
      EXTERNAL MPI_ALLREDUCE
c
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DSQRT

      INTEGER NAEZD

      NAEZD = NAEZ
C     ..
      PI  = 4.D0*ATAN(1.D0)
      FAC = DSQRT((4.0D0*PI)/3.0D0)
      TRP = 0.0D0
      IF (LPOT.LT.1) THEN
         WRITE (6,FMT=9000)
         STOP
 
      END IF
c
C      WRITE (54,FMT=9100)
c
      IREP = 1
c
         IRWS1 = IRWS(IATYP)
         RWS = R(IRWS1,IATYP)
         VOL = 0.25*ALAT**3
c
c---> determine the right potential numbers
c
 
         DO 20 M = -1,1
            LM = 2 + M + 1
c
            DO 30 I = 1,IRWS1
               V1(I) = 0.0d0
   30       CONTINUE
c
            DO 40 ISPIN = 1,NSPIN
c
               IPOT = ISPIN
c
c
               DV = (-3.0D0*V(1,LM,IPOT)-10.0D0*V(2,LM,IPOT)+
     +             18.0D0*V(3,LM,IPOT)-6.0D0*V(4,LM,IPOT)+V(5,LM,IPOT))/
     +              (12.0D0*DRDI(2,IATYP))
c
               V1(2) = RHOC(2,IPOT)* (2.0D0*V(2,LM,IPOT)/R(2,IATYP)+DV)/
     +                 (4.0D0*PI) + V1(2)
c
               DO 50 I = 3,IRWS1 - 2
c
                  DV = (V(I-2,LM,IPOT)-V(I+2,LM,IPOT)+
     +                 8.0D0* (V(I+1,LM,IPOT)-V(I-1,LM,IPOT)))/
     +                 (12.0D0*DRDI(I,IATYP))
c
                  V1(I) = RHOC(I,IPOT)* (2.0D0*V(I,LM,IPOT)/R(I,IATYP)+
     +                    DV)/ (4.0D0*PI) + V1(I)
   50          CONTINUE
c
               DV = (-V(IRWS1-4,LM,IPOT)+6.0D0*V(IRWS1-3,LM,IPOT)-
     +              18.0D0*V(IRWS1-2,LM,IPOT)+10.0D0*V(IRWS1-1,LM,IPOT)+
     +             3.0D0*V(IRWS1,LM,IPOT))/ (12.0D0*DRDI(IRWS1-1,IATYP))
               V1(IRWS1-1) = RHOC(IRWS1-1,IPOT)*
     +                       (2.0D0*V(IRWS1-1,LM,IPOT)/R(IRWS1-1,IATYP)+
     +                       DV)/ (4.0D0*PI) + V1(IRWS1-1)
c
               DV = (3.0D0*V(IRWS1-4,LM,IPOT)-16.0D0*V(IRWS1-3,LM,IPOT)+
     +              36.0D0*V(IRWS1-2,LM,IPOT)-48.0D0*V(IRWS1-1,LM,IPOT)+
     +              25.0D0*V(IRWS1,LM,IPOT))/ (12.0D0*DRDI(IRWS1,IATYP))
c
               V1(IRWS1) = RHOC(IRWS1,IPOT)*
     +                     (2.0D0*V(IRWS1,LM,IPOT)/R(IRWS1,IATYP)+DV)/
     +                     (4.0D0*PI) + V1(IRWS1)
   40       CONTINUE
c
c---> integrate with simpson subroutine
c
            CALL SIMP3(V1,VINT1,1,IRWS1,DRDI(1,IATYP))
c
            FLMH(M,IATYP) = FLM(M,IATYP) - FLMC(M,IATYP)
            FLMXC(M,IATYP) = -FAC*VINT1 - FLMC(M,IATYP)
            FLM(M,IATYP) = FLM(M,IATYP) + FLMXC(M,IATYP)
c
 
   20    CONTINUE
c
C         WRITE (54,FMT=9600) FLMH(1,IATYP),FLMC(1,IATYP),FLMXC(1,IATYP),
C     +     FLM(1,IATYP)
C         WRITE (54,FMT=9601) FLMH(-1,IATYP),FLMC(-1,IATYP),
C     +     FLMXC(-1,IATYP),FLM(-1,IATYP)
C         WRITE (54,FMT=9602) FLMH(0,IATYP),FLMC(0,IATYP),FLMXC(0,IATYP),
C     +     FLM(0,IATYP)
c
        DO I=1,NAEZD
          F(1,I) = 0.0D0
          F(2,I) = 0.0D0
          F(3,I) = 0.0D0
          FALL(1,I) = 0.0D0
          FALL(2,I) = 0.0D0
          FALL(3,I) = 0.0D0
        ENDDO

        F(1,IATYP) = FLM(1,IATYP)
        F(2,IATYP) = FLM(-1,IATYP)
        F(3,IATYP) = FLM(0,IATYP)

C reduce array F to processor with rank 0 to write in correct order to file 'force'

        CALL MPI_ALLREDUCE(F,FALL,3*NAEZD,MPI_DOUBLE_PRECISION,
     +                     MPI_SUM,LCOMM(LMPIC),IERR)

        IF(MYLRANK(LMPIC).EQ.0) THEN
          DO I=1,NAEZD
            WRITE(54,8999) 'force on atom',I,' Z=',ZAT(I),' :',
     +                     FALL(1,I),FALL(2,I),FALL(3,I)
          ENDDO
        ENDIF
c
 8999 FORMAT(A,1X,I5,A,F5.0,A,3(F8.4))

 9000 FORMAT (13x,'error stop in subroutine force :',
     +       ' the charge density has to contain non spherical',
     +       ' contributions up to l=1 at least ')
 9100 FORMAT (1x,33 ('-'),' force on the nucleus ',33 ('-'),/,34x,
     +       ' in units Ry/(a(Bohr) ')
 9600 FORMAT (7x,'fhx=',e12.6,2x,'fcx=',e12.6,2x,'fxcx=',e12.6,2x,'fx=',
     +       e12.6,' Ry/(a(Bohr))')
 9601 FORMAT (7x,'fhy=',e12.6,2x,'fcy=',e12.6,2x,'fxcy=',e12.6,2x,'fy=',
     +       e12.6,' Ry/(a(Bohr))')
 9602 FORMAT (7x,'fhz=',e12.6,2x,'fcz=',e12.6,2x,'fxcz=',e12.6,2x,'fz=',
     +       e12.6,' Ry/(a(Bohr))')
 
      END
