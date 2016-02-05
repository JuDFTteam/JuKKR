c 04.10.95 ***************************************************************
      SUBROUTINE DGDK(GLLKE,DELTAK,DGLLKE,GSPARSE,
     +                 INGSP,NSPBLOCK,ALAT,NAEZ,CLS,EQINV,
     &                 NACLS,RR,EZOA,ATOM,BZKP,IE,KAOEZ,RCLS,GINP)
c
c
c                                      update for sparse matrix 14.6.2001
c   simple to turm of GLLKE array in case of big calculation
c   gsparse can be turned of from inc.p directly when not used
c   set   nauxspd to 1 to turn it off!          
c ************************************************************************
      implicit none
C     .. Parameters ..
      include 'inc.p'
      include 'inc.cls'
      INTEGER LMAX
      PARAMETER (LMAX=LMAXD)
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAX+1)**2)
      INTEGER ALM
      PARAMETER (ALM=LMAXSQ*NDIMGK)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT
      INTEGER IE,NAEZ,NSPBLOCK
C     ..
C     .. Array Arguments ..
      INTEGER 
     +     ATOM(NACLSD,*),
     +     CLS(*),
     +     EZOA(NACLSD,*),
     +     EQINV(*),
     +     KAOEZ(*),
     +     NACLS(*),INGSP(NLSPD,NLSPD)
c
      DOUBLE COMPLEX,intent(in) :: GINP(LMAXSQ*NACLSD,LMAXSQ,NCLSD),
     +                             GLLKE(ALM,ALM),
     &                      GSPARSE(NSPBLOCK*LMAXSQ,NSPBLOCK*LMAXSQ,*)

      DOUBLE COMPLEX ,INTENT(OUT)::  DGLLKE(ALM,ALM)
c
      DOUBLE COMPLEX ::  GLLKEDK(ALM,ALM)

      DOUBLE PRECISION 
     +     BZKP(*),DELTAK(6),
     +     RR(3,0:NRD),
     +     RCLS(3,NACLSD,*)
C     ..
C     .. Local Scalars ..
      INTEGER I,LM1,LM2
c     ..
c     .. Local Arrays ..
      DOUBLE PRECISION KP(6),DKABS,DK_W
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT,DLKE0
C     ..
c ------------------------------------------------------------------------

      CALL CINIT(ALM*ALM,DGLLKE(1,1))
      CALL CINIT(ALM*ALM,GLLKEDK(1,1))

c      WRITE(6,*) "ALM", ALM
c      WRITE(6,*) "LMAXSQ", LMAXSQ

      DKABS=0d0
      DO I=1,3
        DKABS=DKABS+DELTAK(I)**2
      END DO
      DKABS=SQRT(DKABS)
c      DK_W=0.05
      DK_W=DKABS

      DO I=1,3
        KP(I) = BZKP(I)+DELTAK(I)/DKABS*DK_W 
c        KP(I) = BZKP(I)+0.005
      END DO

      DO I=4,6
        KP(I)=0d0
      END DO


      GLLKEDK=0d0

      CALL DLKE0(GLLKEDK,GSPARSE(1:LMAXSQ,1:LMAXSQ,1:NAUXSPD),
     +           INGSP(1:NAEZ,1:NAEZ),1,ALAT,NAEZ,
     +       CLS(1:NAEZ),EQINV(1:NAEZ),NACLS,RR,EZOA(1:NACLSD,1:NAEZ),
     +       ATOM,KP,IE,KAOEZ(1:(NAEZ+NEMBD)),RCLS(:,1:NACLSD,1),
     +       GINP(1:LMAXSQ*NACLSD,1:LMAXSQ,1))

c      WRITE(147,"(3e17.9)") (KP(I),I=1,3)     

      DO LM2=1,ALM
        DO LM1=1,ALM

          DGLLKE(LM1,LM2)=(GLLKEDK(LM1,LM2)-GLLKE(LM1,LM2))
     +                                           /DK_W*DKABS
c         WRITE(147, "((2I5),(6e17.9))") LM2,LM1,GLLKE(LM1,LM2),
c     +                GLLKEDK(LM1,LM2),DGLLKE(LM1,LM2) 
        END DO
      END DO
 
c      WRITE(147, *)

      RETURN

      END                          










