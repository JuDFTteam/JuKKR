      SUBROUTINE CSINWD_SO(F,FINT,DIM_F,IRMIN,IRM,IPAN,IRCUT)
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
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
C      INTEGER IRMIND
c      PARAMETER (IRMIND=IRMD-IRNSD)
      DOUBLE PRECISION A1,A2
      PARAMETER (A1=1.D0/3.D0,A2=4.D0/3.D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER IPAN,IRM,IRMIN,DIM_F
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEXF(DIM_F,DIM_F,IRMIN:IRM),FINT(DIM_F,DIM_F,IRMIN:IRM)
      INTEGER IRCUT(0:IPAN)
C     ..
C     .. Local Scalars ..
      INTEGER I,IEN,IP,IST,LL1,LL2
C     ..
c
c---> loop over kinks
c
      FINT(:,:,:)=0d0
      DO 50 IP = IPAN,1,-1
        IST = IRCUT(IP)
        IEN = IRCUT(IP-1) + 1
        IF (IP.EQ.1) IEN = IRMIN
c        write(6,*) "IST",IST
c        write(6,*) "IEN",IEN
c
        IF (IP.EQ.IPAN) THEN
          DO LL1 = 1,DIM_F
            DO LL2 = 1,DIM_F
            FINT(LL2,LL1,IST) = 0.0D0
c---> integrate fint(ist-1) with a 4 point lagrangian
            FINT(LL2,LL1,IST-1) = (F(LL2,LL1,IST-3) -
     +                       5.0D0*F(LL2,LL1,IST-2)+
     +                      19.0D0*F(LL2,LL1,IST-1)+
     +                       9.0D0*F(LL2,LL1,IST))/24.0D0
            END DO
          END DO
        ELSE
          DO LL1 = 1,DIM_F
            DO LL2 = 1,DIM_F
            FINT(LL2,LL1,IST) = FINT(LL2,LL1,IST+1)
c---> integrate fint(ist-1) with a 4 point lagrangian
            FINT(LL2,LL1,IST-1) = FINT(LL2,LL1,IST+1) +
     +                       (F(LL2,LL1,IST-3)-5.0D0*F(LL2,LL1,IST-2)+
     +                       19.0D0*F(LL2,LL1,IST-1)+
     +                          9.0D0*F(LL2,LL1,IST))/24.0D0
            END DO
          END DO
        END IF
c
c---> calculate fint with an extended 3-point-simpson
c
c        write(123,"((I5),(4e17.9))") IST,F(1,IST)
c        write(123,"((I5),(4e17.9))") IST-1,F(1,IST-1)

        DO 40 I = IST - 2,IEN,-1
          DO LL1 = 1,DIM_F
            DO LL2 = 1,DIM_F
            FINT(LL2,LL1,I) = ((FINT(LL2,LL1,I+2)+F(LL2,LL1,I+2)*A1)+
     +                F(LL2,LL1,I+1)*A2) +
     +                F(LL2,LL1,I)*A1
c          write(123,"((I5),(4e17.9))") I,F(1,I),FINT(1,I)
            END DO
          END DO
   40   CONTINUE
   50 CONTINUE
c
      END
