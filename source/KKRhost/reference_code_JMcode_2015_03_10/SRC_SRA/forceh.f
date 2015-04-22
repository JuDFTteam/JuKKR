      SUBROUTINE FORCEH(CMOM,FLMH,LMAX,NSPIN,NSTART,NEND,R2RHO,V,R,DRDI,
     +                  IRWS,Z)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     calculates the force on nucleus m with hellmann - feynman theorem
c     from a given non spherical charge density at the nucleus site r
c
 
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER LMAX,NEND,NSPIN,NSTART
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CMOM(LMPOTD,*),DRDI(IRMD,*),FLMH(-1:1,*),
     +       R(IRMD,*),
     +       R2RHO(IRMD,LMPOTD,NATYPD,*),V(IRMD,LMPOTD,*),Z(*)
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
      IF (LMAX.LT.1) THEN
         WRITE (6,FMT=9000)
         STOP
 
      END IF
c
c---> loop over the rep. atoms
c
      DO 10 IATYP = NSTART,NEND
c
c---> reading the right Wigner-S. radius
c
         IRWS1 = IRWS(IATYP)
         RWS = R(IRWS1,IATYP)
c
c---> determine the right potential numbers
c
         IPOT = NSPIN* (IATYP-1) + 1
 
         DO 20 M = -1,1
            LM = 2 + M + 1
c
            V1(1) = 0.0D0
            DO 30 I = 2,IRWS1
               V1(I) = R2RHO(I,LM,IATYP,1)* (R(I,IATYP)** (-2.0D0))
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
     +                 2.0D0*CMOM(LM,IATYP)/ (RWS**3)
c
c---> total Hellman-Feynman force
c
            FLMH(M,IATYP) = (FLM(M,1)+FLM(M,2))*Z(IATYP)
   20    CONTINUE
   10 CONTINUE
c
c
 9000 FORMAT (13x,'error stop in subroutine force :',
     +       ' the charge density has to contain non spherical',
     +       ' contributions up to l=1 at least ')
 
      END
