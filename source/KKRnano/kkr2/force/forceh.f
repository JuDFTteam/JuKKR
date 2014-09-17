C     Calculates the Hellman-Feynman force.
C
C     ATTENTION: a factor sqrt(4*pi/3) is missing here, this is
C     corrected in routine 'FORCE'
      SUBROUTINE FORCEH(FLMH,LPOT,RHO2NS,V,R,DRDI,
     &                  IRWS,Z,irmd)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     calculates the force on nucleus m with hellmann - feynman theorem
c     from a given non spherical charge density at the nucleus site r
c
 
c-----------------------------------------------------------------------
C     .. Parameters ..

      INTEGER irmd

C     .. Scalar Arguments ..
      INTEGER LPOT

      DOUBLE PRECISION DRDI(IRMD),FLMH(-1:1),
     &       R(IRMD),
     &       RHO2NS(IRMD,(LPOT+1)**2), V(IRMD,(LPOT+1)**2,2), Z
      INTEGER IRWS
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION PI,RWS,VINT1
      INTEGER I,IPOT,IRWS1,LM,M
      integer ind
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

      IRWS1 = IRWS
      RWS = R(IRWS1)

      do ind = 1, IRWS

         RWS = R(ind)
         IRWS1 = ind

         IPOT = 1
 
         DO 20 M = -1,1
            LM = 2 + M + 1
c
            V1(1) = 0.0D0
            DO 30 I = 2,IRWS1
               V1(I) = RHO2NS(I,LM)* (R(I)** (-2.0D0) - R(I) / RWS**3)
   30       CONTINUE
c
c---> integrate with simpson subroutine
c
            CALL SIMP3(V1,VINT1,1,IRWS1,DRDI)
c
            FLM(M,1) = 2.0D0*VINT1
c
c---> use coulomb potential to determine extra atomic contribution
c
            FLM(M,2) = V(IRWS1,LM,IPOT)* (3.0D0/ (4.0D0*PI*RWS))
c
c---> total Hellman-Feynman force
c
            FLMH(M) = (FLM(M,1)+FLM(M,2))*Z
   20    CONTINUE


      write(*, '(2(e23.17, X))') R(ind),
     &     FLMH(1) * sqrt(16.0d0 * atan(1.0d0)/3.0d0)
      enddo
c
c
 9000 FORMAT (13x,'error stop in subroutine force :',
     +       ' the charge density has to contain non spherical',
     +       ' contributions up to l=1 at least ')
 
      END
