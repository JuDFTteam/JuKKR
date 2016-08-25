      SUBROUTINE CSIMPK1(CF,CFINT,LMMSQD,IRMIN,IRM,DRDI)
c ************************************************************************
c     this subroutine does an integration up to rcut of an
c     complex function cf with an extended 3-point-simpson :
c
c                                  rcut
c                      cfint(ll) = { cf(ll,r') dr'
c                                  0
c
c     modified for functions with kinks - at each kink the
c     integration is restarted .
c
c     attention : input cf is destroyed !
c
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
C     ..
C     .. Scalar Arguments ..
      INTEGER IRMIN,IRM,LMMSQD
      DOUBLE COMPLEX CFINT(LMMSQD,IRM)
c      DOUBLE COMPLEX CFINT
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX CF(LMMSQD,IRM)
      DOUBLE PRECISION DRDI(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,A2
      INTEGER I,IEN,IP,IST,N,LL
C     ..
C     .. External Functions ..
      DOUBLE COMPLEX CSUM
      EXTERNAL CSUM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
      A1 = 1.0D0/3.0D0
      A2 = 4.0D0/3.0D0
c      CFINT = 0.0D0
c
c
c---> loop over kinks
c
        IST = IRMIN
        IEN = IRM
c        IST = IRCUT(IP-1)+1
c        IEN = IRCUT(IP)
c
        DO 10 I = IST,IEN
          DO LL=1,LMMSQD
           CF(LL,I) = CF(LL,I)*DRDI(I)
          ENDDO
   10   CONTINUE
c
c        IF (MOD(IEN-IST,2).EQ.0) THEN
          DO LL=1,LMMSQD
          CFINT(LL,IST) = 0d0
          CFINT(LL,IST+1)=(CF(LL,IST+3)-5D0*CF(LL,IST+2)+
     +                    19D0*CF(LL,IST+1)+9D0*CF(LL,IST))/24D0
          ENDDO
c          IST = IST + 1
c          N = (IEN-IST+1)/2

c        ELSE
c---> four point lagrange integration for the first step
c          DO LL=1,LMMSQD
c          CFINT(LL) = CFINT(LL)+(9.0D0*CF(LL,IST)+
c     +            19.0D0*CF(LL,IST+1)-
c     +            5.0D0*CF(LL,IST+2)+CF(LL,IST+3))/24.0D0 +
c     +            (CF(LL,IST+1)-CF(LL,IEN))/3.0D0
c          ENDDO
c          IST = IST + 2
c          N = (IEN-IST+1)/2
c        END IF
c
c---> calculate with an extended 3-point-simpson
c
        DO I=IST+2,IEN
        DO LL=1,LMMSQD
        CFINT(LL,I) = CFINT(LL,I-2)+A1*CF(LL,I-2)+
     +                   A2*CF(LL,I-1)+A1*CF(LL,I)
c        CFINT(LL) = CFINT(LL)+A1*CSUM(N,CF(LL,IST),2)+
c     +                   A2*CSUM(N,CF(LL,IST+1),2)
        ENDDO
        ENDDO
c
      END
