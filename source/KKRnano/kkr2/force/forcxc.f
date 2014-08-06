      SUBROUTINE FORCXC(FLM,FLMC,LPOT,NSPIN,RHOC,V,R,
     +                  DRDI,IRWS, irmd)
C
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     calculates the force on nucleus m
c     from a given non spherical charge density at the nucleus site r
c     with core correction(exchange contribution)
 
c-----------------------------------------------------------------------

      INTEGER irmd
C
C      INTEGER LMPOTD
C      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER LPOT,NSPIN
C     ..
C     .. Array Arguments ..

      DOUBLE PRECISION DRDI(IRMD),FLM(-1:1),FLMC(-1:1),R(IRMD),
     &       RHOC(IRMD,*),
     &       V(IRMD,(LPOT+1)**2,2)

      INTEGER IRWS
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DV,FAC,PI,RWS,VINT1
      INTEGER I,IPOT,IREP,IRWS1,ISPIN,LM,M
C     ..
C     .. Local Arrays ..

      DOUBLE PRECISION FLMXC(-1:1),V1(IRMD)
      DOUBLE PRECISION TAIL_COR(-1:1)
C     ..
C     .. External Subroutines ..
      EXTERNAL SIMP3
C     ..
C     .. Save statement ..
      SAVE PI

C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DSQRT
C     ..
      PI  = 4.D0*ATAN(1.D0)
      FAC = DSQRT((4.0D0*PI)/3.0D0)

      IF (LPOT.LT.1) THEN
         WRITE (6,FMT=9000)
         STOP
 
      END IF

         TAIL_COR = 0.0d0

         IRWS1 = IRWS
         RWS = R(IRWS1)
 
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
     +              (12.0D0*DRDI(2))
c
               V1(2) = RHOC(2,IPOT)* (2.0D0*V(2,LM,IPOT)/R(2)+DV)/
     +                 (4.0D0*PI) + V1(2)
c
               DO 50 I = 3,IRWS1 - 2
c
                  DV = (V(I-2,LM,IPOT)-V(I+2,LM,IPOT)+
     +                 8.0D0* (V(I+1,LM,IPOT)-V(I-1,LM,IPOT)))/
     +                 (12.0D0*DRDI(I))
c
                  V1(I) = RHOC(I,IPOT)* (2.0D0*V(I,LM,IPOT)/R(I)+
     +                    DV)/ (4.0D0*PI) + V1(I)
   50          CONTINUE
c
               DV = (-V(IRWS1-4,LM,IPOT)+6.0D0*V(IRWS1-3,LM,IPOT)-
     +              18.0D0*V(IRWS1-2,LM,IPOT)+10.0D0*V(IRWS1-1,LM,IPOT)+
     +             3.0D0*V(IRWS1,LM,IPOT))/ (12.0D0*DRDI(IRWS1-1))
               V1(IRWS1-1) = RHOC(IRWS1-1,IPOT)*
     +                       (2.0D0*V(IRWS1-1,LM,IPOT)/R(IRWS1-1)+
     +                       DV)/ (4.0D0*PI) + V1(IRWS1-1)
c
               DV = (3.0D0*V(IRWS1-4,LM,IPOT)-16.0D0*V(IRWS1-3,LM,IPOT)+
     +              36.0D0*V(IRWS1-2,LM,IPOT)-48.0D0*V(IRWS1-1,LM,IPOT)+
     +              25.0D0*V(IRWS1,LM,IPOT))/ (12.0D0*DRDI(IRWS1))
c
               V1(IRWS1) = RHOC(IRWS1,IPOT)*
     +                     (2.0D0*V(IRWS1,LM,IPOT)/R(IRWS1)+DV)/
     +                     (4.0D0*PI) + V1(IRWS1)

CDEBUG      Tail correction
            TAIL_COR(M) = TAIL_COR(M) +
     &                    FAC * RHOC(IRWS1, IPOT) / (4.0d0*PI) *
     &                    V(IRWS1, LM, IPOT)
CDEBUG

   40       CONTINUE
c
c---> integrate with simpson subroutine
c
            CALL SIMP3(V1,VINT1,1,IRWS1,DRDI)

            FLMXC(M) = -FAC*VINT1 - FLMC(M)
            FLM(M) = FLM(M) + FLMXC(M)

   20    CONTINUE

CDEBUG Tail correction
C      FLM = FLM + TAIL_COR
CDEBUG

C     Result: total force in FLM

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
