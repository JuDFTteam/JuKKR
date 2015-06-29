      SUBROUTINE CSINWD(F,FINT,LMMSQD,IRMIND,IRMD,IPAN,IRCUT)
c-----------------------------------------------------------------------
c     this subroutine does an inwards integration of llmax
c     functions f with an extended 3-point-simpson :
c
c
c                               irmax
c                   fint(ll,i) = { f(ll,i') di'
c                                ir
c
c  the starting value for this integration at ist - 1 is determined by
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
c    modified by m. ogura, june 2015
c-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION A1,A2,A3
      PARAMETER (A1=5.D0/12.D0,A2=8.D0/12.D0,A3=-1.D0/12.D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER IPAN,IRMD,IRMIND,LMMSQD
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX F(LMMSQD,IRMIND:IRMD),FINT(LMMSQD,IRMIND:IRMD)
      INTEGER IRCUT(0:IPAN)
C     ..
C     .. Local Scalars ..
      INTEGER I,IEN,IP,IST,LL
C     ..
c
c---> loop over kinks
c
      DO 50 IP = IPAN,1,-1
        IST = IRCUT(IP)
        IEN = IRCUT(IP-1) + 1
        IF (IP.EQ.1) IEN = IRMIND
c
        IF (IP.EQ.IPAN) THEN
          DO 10 LL = 1,LMMSQD
            FINT(LL,IST) = 0.0D0
   10     CONTINUE

        ELSE
          DO 20 LL = 1,LMMSQD
            FINT(LL,IST) = FINT(LL,IST+1)
   20     CONTINUE
        END IF
c
c---> calculate fint with an extended 3-point-simpson
c
        DO 40 I = IST,IEN+2,-2
           DO 30 LL = 1,LMMSQD
              FINT(LL,I-1)=FINT(LL,I)
     +                    +F(LL,I)*A1+F(LL,I-1)*A2+F(LL,I-2)*A3
              FINT(LL,I-2)=FINT(LL,I-1)
     +                    +F(LL,I)*A3+F(LL,I-1)*A2+F(LL,I-2)*A1
 30        CONTINUE
 40     CONTINUE
        IF(MOD(IST-IEN,2).EQ.1)THEN
           DO 60 LL=1,LMMSQD
 60           FINT(LL,IEN)=FINT(LL,IEN+1)
     +                    +F(LL,IEN)*A1+F(LL,IEN+1)*A2+F(LL,IEN+2)*A3
           ENDIF
 50     CONTINUE
c     
      END
