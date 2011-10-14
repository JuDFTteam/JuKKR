      SUBROUTINE FORCEH(CMOM,FLMH,LPOT,NSPIN,IATYP,RHO2NS,V,R,DRDI,
     &                  IRWS,Z,
C                       new input parameter after inc.p replace
     &                  irmd)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     calculates the force on nucleus m with hellmann - feynman theorem
c     from a given non spherical charge density at the nucleus site r
c
 
c-----------------------------------------------------------------------
C     .. Parameters ..

      INTEGER irmd

C     INTEGER LMPOTD
C     PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER LPOT,NSPIN
C     ..
C     .. Array Arguments ..
C     DOUBLE PRECISION CMOM(LMPOTD),DRDI(IRMD,*),FLMH(-1:1,*),
C    +       R(IRMD,*),
C    +       RHO2NS(IRMD,LMPOTD),V(IRMD,LMPOTD,2),Z(*)

      DOUBLE PRECISION CMOM((LPOT+1)**2),DRDI(IRMD,*),FLMH(-1:1,*),
     &       R(IRMD,*),
     &       RHO2NS(IRMD,(LPOT+1)**2), V(IRMD,(LPOT+1)**2,2), Z(*)
      INTEGER IRWS(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION PI,RWS,VINT1
      INTEGER I,IATYP,IPOT,IRWS1,LM,M
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION FLM(-1:1,2),V1(IRMD)
C     ..
C     .. External Subroutines ..
      EXTERNAL SIMP3
C     ..
C     .. Save statement ..
      SAVE PI
C     ..
c
C     .. Intrinsic Functions ..
      INTRINSIC ATAN
C     ..
      PI = 4.D0*ATAN(1.D0)
      IF (LPOT.LT.1) THEN
         WRITE (6,FMT=9000)
         STOP
 
      END IF
c
c
c---> reading the right Wigner-S. radius
c
         IRWS1 = IRWS(IATYP)
         RWS = R(IRWS1,IATYP)
c
c---> determine the right potential numbers
c
         IPOT = 1
 
         DO 20 M = -1,1
            LM = 2 + M + 1
c
            V1(1) = 0.0D0
            DO 30 I = 2,IRWS1
               V1(I) = RHO2NS(I,LM)* (R(I,IATYP)** (-2.0D0))
   30       CONTINUE
c
c---> integrate with simpson subroutine
c
            CALL SIMP3(V1,VINT1,1,IRWS1,DRDI(1,IATYP))
c
            FLM(M,1) = 2.0D0*VINT1
c
c---> use coulomb potential to determine extra atomic contribution
c
            FLM(M,2) = V(IRWS1,LM,IPOT)* (3.0D0/ (4.0D0*PI*RWS)) -
     +                 2.0D0*CMOM(LM)/ (RWS**3)
c
c---> total Hellman-Feynman force
c
            FLMH(M,IATYP) = (FLM(M,1)+FLM(M,2))*Z(IATYP)
   20    CONTINUE
c
c
 9000 FORMAT (13x,'error stop in subroutine force :',
     +       ' the charge density has to contain non spherical',
     +       ' contributions up to l=1 at least ')
 
      END
