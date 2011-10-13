c 14.10.95 *************************************************************
      SUBROUTINE SOUTK(F,FINT,IPAN,IRCUT)
c ************************************************************************
c    this subroutine does an outwards integration of a function
c    with kinks
c
c
c                             r
c                   fint(r) = s f(r') dr'
c                             0
c
c    at each kink the integration is restarted
c    the starting value for this integration is determined by
c    a 4 point lagrangian integration  , coefficients given by
c    m. abramowitz and i.a. stegun, handbook of mathematical functions,
c    nbs applied mathematics series 55 (1968)
c
c    the weights drdi have to be multiplied before calling this
c    subroutine .
c
c                                     b. drittler oct. 1989
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
C     ..
C     .. Scalar Arguments ..
      INTEGER IPAN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION F(*),FINT(*)
      INTEGER IRCUT(0:IPAND)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,A2
      INTEGER I,IEN,IP,IST
C     ..
      A1 = 1.0D0/3.0D0
      A2 = 4.0D0/3.0D0
c
c---> loop over kinks
c
      DO 20 IP = 1,IPAN
        IEN = IRCUT(IP)
        IST = IRCUT(IP-1) + 1
c
        IF (IP.EQ.1) THEN
          FINT(IST) = 0.0D0
c---> integrate fint(ist+1) with a 4 point lagrangian
          FINT(IST+1) = (F(IST+3)-5.0D0*F(IST+2)+19.0D0*F(IST+1)+
     +                  9.0D0*F(IST))/24.0D0

        ELSE
          FINT(IST) = FINT(IST-1)
c---> integrate fint(ist+1) with a 4 point lagrangian
          FINT(IST+1) = FINT(IST-1) + (F(IST+3)-5.0D0*F(IST+2)+
     +                  19.0D0*F(IST+1)+9.0D0*F(IST))/24.0D0
        END IF

c
c---> calculate fint with an extended 3-point-simpson
c
        DO 10 I = IST + 2,IEN
          FINT(I) = ((FINT(I-2)+F(I-2)*A1)+F(I-1)*A2) + F(I)*A1
   10   CONTINUE
   20 CONTINUE
c
      END
