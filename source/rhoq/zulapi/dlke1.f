c 01.06.99 ***************************************************************
      SUBROUTINE DLKE1(GLLKE,ALAT,NACLS,RR,EZOA,
     +                 ATOM,BZKP,IE,IC,FAC,GINP,RCLS)
      implicit none
c ************************************************************************
c
c     Fourier transformation of the cluster Greens function GINP
c
c ------------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
      include 'inc.cls'
c      INTEGER NATOMD
c      PARAMETER (NATOMD=79)
c      INTEGER LMAX
c      PARAMETER (LMAX=4)
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAXD+1)**2)
      INTEGER ALM,CLM
      PARAMETER (ALM=LMAXSQ*NAEZD,CLM=LMAXSQ*NACLSD)
      DOUBLE COMPLEX CZERO,CI
      PARAMETER (CZERO= (0.0D0,0.0D0),CI= (0.0D0,1.0D0))
      DOUBLE COMPLEX CONE,CONEM
      PARAMETER (CONE=(1.0D0,0.0D0),CONEM=(-1.0D0,0.0D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT
      INTEGER IC,IE
      DOUBLE COMPLEX FAC
C     ..
C     .. Array Arguments ..
      INTEGER ATOM(*),
     +        EZOA(*),
     +        NACLS
      DOUBLE COMPLEX, intent(out):: GLLKE(ALM,LMAXSQ)
      DOUBLE COMPLEX, intent(in) :: GINP(LMAXSQ*NACLSD,LMAXSQ)
      DOUBLE PRECISION BZKP(3),
     +                 RR(3,0:NRD),
     +                 RCLS(3,*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CONVPU,R,TPI
      INTEGER AM,I,II,IESAVE,IM,LM2,M,ML,N,N1,NAC,NL,J
      DOUBLE COMPLEX  EIKR,FACNM,TT
      LOGICAL LSTART,OPT,TEST
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX ARG(3),
     +               GMN(LMAXSQ*NACLSD,LMAXSQ),
     &               CONEMN(LMAXSQ,LMAXSQ)
      INTEGER LF(144)
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT,TEST,OPT,ZAXPY
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DATAN,EXP
C     ..
C     .. COMMON BLOCK
      DOUBLE PRECISION RFOURIER
      COMMON /RFOUR/ RFOURIER
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Data statements ..
      DATA LF/0,3*1,5*2,7*3,9*4,11*5,13*6,15*7,17*8,19*9,21*10,23*11/
c      DATA LF / 0,3*1,5*2,7*3,9*4,11*5,13*6 /
      DATA LSTART /.TRUE./
c ------------------------------------------------------------------------
C
c      WRITE(6,*) "In dlke1"

      GLLKE=0d0
      GMN=0d0
      II = 3
      IF (OPT('COMPLEX ')) II = 6
      IF (TEST('BZKP    ')) write(6,FMT='(6f12.6)') (bzkp(i),i=1,ii)
c
      IF (LSTART) THEN
c        IESAVE = 0
c        OPEN (17,FILE=I17,FORM='unformatted')
c        READ (69) ((RM(I,M),I=1,3),M=1,NATOM)
c        WRITE (6,FMT=*) ((RM(I,M),I=1,3),M=1,NATOM)
      DO M = 1,LMAXSQ
        DO N = 1,LMAXSQ
        CONEMN(M,N)=CONEM**(LF(M)+LF(N))
        END DO
      END DO
        LSTART = .FALSE.
      END IF

      CONVPU = ALAT/2.D0/3.14159265358979312D0
      TPI = 8.0D0*DATAN(1.0D0)                 ! = 2*PI

      DO 40 M = 1,LMAXSQ
        DO 30 N = 1,LMAXSQ
          FACNM = FAC**(LF(N)+LF(M))
          DO 20 N1 = 1,NACLSD
            NL = (N1-1)*LMAXSQ + N
            ML = (N1-1)*LMAXSQ + M
c            GMN(NL,M) = FACNM*(0.5D0,0.0d0)*
c     +           ( CONEMN(M,N)*GINP(ML,N)+GINP(NL,M) )
c  Changed 26.11.2001 

C            GMN(NL,M) = (CONEM**(LF(M)+LF(N))*GINP(ML,N)+
c     +           GINP(NL,M))*0.5D0

cccc correct
            GMN(NL,M) = FACNM*GINP(NL,M)
c            WRITE(67,"((2I5),(2e17.9))") M,NL,GMN(NL,M)
cccc correct
 20       CONTINUE
 30     CONTINUE
 40   CONTINUE

c      DO 60 M = 1,LMAXSQ
c        DO 50 N = 1,ALM
c          GLLKE(N,M) = CZERO
c 50     CONTINUE
c 60   CONTINUE
      
      CALL CINIT(LMAXSQ*NAEZD*LMAXSQ,GLLKE)

      DO 90 M = 1,NACLS
c      DO 90 M = 1,NACLS(IC)
c        WRITE(166,"(I5,3e17.9)") M,(RCLS(J,M),J=1,3)
c        write(16,"(2I5,3e17.9)") m,ezoa(m),(RR(J,EZOA(M)),J=1,3)

        R = DSQRT(RCLS(1,M)*RCLS(1,M) +
     +            RCLS(2,M)*RCLS(2,M) +
     +            RCLS(3,M)*RCLS(3,M)   )
c
c --->  test for restricted Fourier transformation
c
        IF (RFOURIER.GT.1D-6 .AND. R.GT.RFOURIER+1.d-3) GOTO 90
c
c --->  for option 'WIRE': avoid artifical couplings in the
c       structural Greens Function in in-plane-direction (perp. to c-axis)
c
        IF (ATOM(M).LT.0) GOTO 90
c
        IF (OPT('ONEBULK ')) THEN     ! added 1.02.2000
c                                       corrected on 25.02.2000
c     if the phase factor exp(ik(r-r')) is included      ~
c     in the G...so if we resolve the dyson part for the G
c     and not for the G (see Peter Lang Ph.D thesis)
c     
c     Here we do   --                           nn'
c                  \                            ii'          ii'
c                  /  exp( ik(X  -X + R  -R  ))G   (E)  =   G   (k,E)
c                  --          n'  n   i'  i    LL'          LL'
c                  n'
c                   
c In this case rcls is always (by constraction symmetric around each
c atom this means that a minus sign will not affect the result of the
c summation

           ARG(1) = -CI*TPI*RCLS(1,M)
           ARG(2) = -CI*TPI*RCLS(2,M)
           ARG(3) = -CI*TPI*RCLS(3,M)
        ELSE
c     
c     Here we do   --                  nn'
c                  \                   ii'          ii'
c                  /  exp(+ik(x  -x ))G   (E)  =   G   (k,E)
c                  --          n'  n   LL'          LL'
c                  n'
c  Be carefull a minus sign must be included here. RR is not
c  symmetric around each atom. The minus comes from the fact that
c  the repulsive potential GF is calculated for 0n and not n0!          
c  and that is why we nead a minus sign extra!
c  
           ARG(1) = -CI*TPI*RR(1,EZOA(M))
           ARG(2) = -CI*TPI*RR(2,EZOA(M))
           ARG(3) = -CI*TPI*RR(3,EZOA(M))
c        write(6,888) m,ezoa(m),RR(1,EZOA(M)),RR(2,EZOA(M)),
c     &   RR(3,EZOA(M))  
c 888    format(2I7,3F10.5)
c            write(6,*) 'In fourier ', m,EZOA(M)
c
c Added 8.12.2001  test
c
c           ARG(1) = -CI*TPI*RCLS(1,M)
c           ARG(2) = -CI*TPI*RCLS(2,M)
c           ARG(3) = -CI*TPI*RCLS(3,M)
cc
c  Added 8.12.2001  test
c

        END IF
c
c        write(17,"(I5,6e17.9)") m,(ARG(J),J=1,3)
c        write(18,"(6e17.9)") (BZKP(J),J=1,3)
        TT = BZKP(1)*ARG(1)+BZKP(2)*ARG(2)+BZKP(3)*ARG(3)
c        write(35,*) M,(BZKP(J),J=1,2)
c        write(36,*) M,BZKP(1),TT
c
        IF (OPT('COMPLEX ')) THEN
c          TT = TT + CI*(BZKP(4)*ARG(1)+BZKP(5)*ARG(2)+BZKP(6)*ARG(3))
        END IF
c
c        write(6,*) BZKP(1),ARG(1),BZKP(2),ARG(2),BZKP(3),ARG(3)
c        write(6,*) BZKP(4),BZKP(5),BZKP(6)
c        write(6,*) 'm,atom(m),tt',m,atom(m),tt
c
        EIKR = EXP(TT) * CONVPU    ! convert to p.u.

c        write(6,*) 'eikr',eikr

        IM = 1 + (M-1)      *LMAXSQ
        AM = 1 + (ATOM(M)-1)*LMAXSQ
        DO 80 LM2 = 1,LMAXSQ
          CALL ZAXPY(LMAXSQ,EIKR,GMN(IM,LM2),1,GLLKE(AM,LM2),1)  
   80   CONTINUE

 90   CONTINUE                      ! M = 1,NACLS(IC)

      RETURN
 9000 format(3f12.4)
 9010 format(2f18.10)
      END




