      SUBROUTINE CSOUT_SPLINE(F,FINT,LMMSQD,IRMIN,IRM,IPAN,IRCUT)
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
      use mod_splint, only: splint_complex
      implicit none
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
      INTEGER I,IEN,IP,IST,LL,
     +        NMAX
C     ..
c      ..Local
c      DOUBLE COMPLEX F_NEW(LMMSQD,IRMIN:IRM),FINT_NEW(LMMSQD,IRMIN:IRM)
      DOUBLE COMPLEX, allocatable:: FINT_SP(:,:),FDERSP(:)!,FINTSP(:),
c     +                              FDER_SP(:,:),FEND(:),
      DOUBLE COMPLEX           ::   FDUM
c 
      REAL,allocatable  ::      I_SP(:)  
c---> loop over kinks
c

      IF ( MOD((IRCUT(1)-IRCUT(0)-1),2) .EQ. 1) 
     +                       STOP "Number of points not even"    
      NMAX=(IRCUT(1)-IRMIN)/2+1

c      ALLOCATE(FEND(LMMSQD))
      ALLOCATE(FINT_SP(NMAX,LMMSQD))
c      ALLOCATE(FINT_SP(LMMSQD,NMAX))
      ALLOCATE(FDERSP(NMAX))
c        ALLOCATE(FINTSP(NMAX))
c      ALLOCATE(FDER_SP(LMMSQD,NMAX))
      ALLOCATE(I_SP(NMAX))

      DO 50 IP = 1,IPAN
        IEN = IRCUT(IP)
        IST = IRCUT(IP-1) + 1
c
        IF (IP.EQ.1) THEN
c          IST = 276
          IST = IRMIN
          DO 10 LL = 1,LMMSQD
            FINT(LL,IST) = 0.0D0
c---> integrate fint(ist+1) with a 4 point lagrangian
            FINT(LL,IST+1) = (F(LL,IST+3)-5.0D0*F(LL,IST+2)+
     +                       19.0D0*F(LL,IST+1)+9.0D0*F(LL,IST))/24.0D0
c            FINT(LL,IST+1) = 0.5d0*(F(LL,IST+1)-F(LL,IST))
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

c     implemented by Swantje
c---> for the first panel:                                     
c     interpolate every second point of the integration with a 
c     spline interpolation

        IF (IP.EQ.1) THEN

          DO I = IST,IEN,2
            I_SP((I-IST)/2+1)=I
            DO LL = 1,LMMSQD
              FINT_SP((I-IST)/2+1,LL) = FINT(LL,I)
            END DO   
          END DO  
          DO LL = 1,LMMSQD
                   
            CALL SPLINE(NMAX,I_SP,FINT_SP(1:NMAX,LL),NMAX,F(LL,IST),
     +                            F(LL,IEN),FDERSP)
            DO I = IST+1,IEN-1,2
               
              CALL splint_complex(I_SP,FINT_SP(1:NMAX,LL),FDERSP,NMAX,
     +                             REAL(I),FINT(LL,I),FDUM)
            END DO
          END DO

        END IF

   50 CONTINUE

      DEALLOCATE(FINT_SP)
      DEALLOCATE(FDERSP)
      DEALLOCATE(I_SP)

      END
