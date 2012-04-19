c 14.10.95 *************************************************************
      SUBROUTINE SIMPK(F,FINT,IPAN,IRCUT,DRDI)
c **********************************************************************
c     this subroutine does an integration up to rcut of an
c     real function f with an extended 3-point-simpson :
c
c                            rcut
c                      fint = { f(r') dr'
c                             0
c
c     modified for functions with kinks - at each kink the
c     integration is restarted .
c
c     attention : input f is destroyed !
c
c-----------------------------------------------------------------------

      IMPLICIT NONE

C     .. Scalar Arguments ..
      DOUBLE PRECISION FINT
      INTEGER IPAN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(*),F(*)
      INTEGER IRCUT(0:IPAN)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,A2
      INTEGER I,IEN,IP,IST,N
C     ..
C     .. External Functions ..
      DOUBLE PRECISION SSUM
      EXTERNAL SSUM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
      A1 = 4.0D0/3.0D0
      A2 = 2.0D0/3.0D0
      FINT = 0.0D0
c
      DO 20 IP = 1,IPAN
c
c---> loop over kinks
c
        IST = IRCUT(IP-1) + 1
        IEN = IRCUT(IP)
c
        DO 10 I = IST,IEN
          F(I) = F(I)*DRDI(I)
   10   CONTINUE
        IF (MOD(IEN-IST,2).EQ.0) THEN
          FINT = FINT + (F(IST)-F(IEN))/3.0D0
          IST = IST + 1
          N = (IEN-IST+1)/2

        ELSE
c---> four point lagrange integration for the first step
          FINT = FINT + (9.0D0*F(IST)+19.0D0*F(IST+1)-5.0D0*F(IST+2)+
     +           F(IST+3))/24.0D0 + (F(IST+1)-F(IEN))/3.0D0
          IST = IST + 2
          N = (IEN-IST+1)/2
        END IF
c
c---> calculate with an extended 3-point-simpson
c
        FINT = FINT + A1*SSUM(N,F(IST),2) + A2*SSUM(N,F(IST+1),2)
   20 CONTINUE
c
      END
