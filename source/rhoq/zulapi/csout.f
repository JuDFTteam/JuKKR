      SUBROUTINE CSOUT(F,FINT,LMMSQD,IRMIN,IRM,IPAN,IRCUT)
c-----------------------------------------------------------------------
c     this subroutine does an outwards integration of llmax
c     functions f with an extended 3-point-simpson :
c
c
c                                ir
c                   fint(ll,i) = { f(ll,i') di'
c                               irmin
c
c  the starting value for this integration at irmin+1 is determined by
c    a 4 point lagrangian integration  , coefficients given by
c    m. abramowitz and i.a. stegun, handbook of mathematical functions,
c    nbs applied mathematics series 55 (1968)
c
c  attention in case of radial integration :
c       the weights drdi have to be multiplied before calling this
c       subroutine .
c
c                                     b. drittler mar. 1989
c
c    modified for functions with kinks - at each kink the integration
c      is restarted
c
c    attention : it is supposed that irmin + 3 is less than imt !
c
c
c                                     b. drittler july 1989
c-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION A1,A2
      PARAMETER (A1=1.D0/3.D0,A2=4.D0/3.D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER IPAN,IRM,IRMIN,LMMSQD
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX F(LMMSQD,IRMIN:IRM),FINT(LMMSQD,IRMIN:IRM)
      INTEGER IRCUT(0:IPAN)
C     ..
C     .. Local Scalars ..
      INTEGER I,IEN,IP,IST,LL
C     ..
c
c---> loop over kinks
c
      DO 50 IP = 1,IPAN
        IEN = IRCUT(IP)
        IST = IRCUT(IP-1) + 1
c
        IF (IP.EQ.1) THEN
          IST = IRMIN
          DO 10 LL = 1,LMMSQD
            FINT(LL,IST) = 0.0D0
c---> integrate fint(ist+1) with a 4 point lagrangian
            FINT(LL,IST+1) = (F(LL,IST+3)-5.0D0*F(LL,IST+2)+
     +                       19.0D0*F(LL,IST+1)+9.0D0*F(LL,IST))/24.0D0
   10     CONTINUE

        ELSE
          DO 20 LL = 1,LMMSQD
            FINT(LL,IST) = FINT(LL,IST-1)
c---> integrate fint(ist+1) with a 4 point lagrangian
            FINT(LL,IST+1) = FINT(LL,IST-1) +
     +                       (F(LL,IST+3)-5.0D0*F(LL,IST+2)+
     +                       19.0D0*F(LL,IST+1)+9.0D0*F(LL,IST))/24.0D0
   20     CONTINUE
        END IF

c
c---> calculate fint with an extended 3-point-simpson
c
        DO 40 I = IST + 2,IEN
          DO 30 LL = 1,LMMSQD
            FINT(LL,I) = ((FINT(LL,I-2)+F(LL,I-2)*A1)+F(LL,I-1)*A2) +
     +                   F(LL,I)*A1
   30     CONTINUE
   40   CONTINUE
   50 CONTINUE
c
      END
