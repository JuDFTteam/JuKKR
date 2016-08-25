c ************************************************************************
      SUBROUTINE TAUTOG1_SO(GS,TINVLL,DSYMLL,NSHELL,RFCTOR,GLL0,IGF,
     +                   TAUVBZ,NSYMAT,NSH1,NSH2,RATOM)
      implicit none
c ************************************************************************
c
c     GLL0 = GMATLL in subroutine GLL96
c
c ------------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
c
      INTEGER NSYMAXD
      PARAMETER (NSYMAXD=48)
      DOUBLE COMPLEX CZERO,CONE
      PARAMETER (CZERO= (0.0D0,0.0D0),CONE= (1.D0,0.D0))
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAXD+1)**2)
      INTEGER LMMAXSO
      PARAMETER (LMMAXSO= NSPOD*LMAXSQ)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX TAUVBZ
      DOUBLE PRECISION RFCTOR
      INTEGER IGF,NSYMAT,NSHELL
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX 
     +     GLL0(LMMAXSO,LMMAXSO,*),
     +     GS(LMMAXSO,LMMAXSO,NSYMAT,*),
     +     TINVLL(LMMAXSO,LMMAXSO,*),
     +     DSYMLL(LMMAXSO,LMMAXSO,*)
      DOUBLE PRECISION RATOM(3,*)
      INTEGER NSH1(*),NSH2(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,IND,IU,LM1,LM2,NSLL,NS,J
      LOGICAL LDIA
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX GLL(LMMAXSO,LMMAXSO),
     +               TPG(LMMAXSO,LMMAXSO),
     +               XC(LMMAXSO,LMMAXSO)
c      DOUBLE COMPLEX GLL(LMAXSQ,LMAXSQ),
c     +               TPG(LMAXSQ,LMAXSQ),
c     +               XC(LMAXSQ,LMAXSQ)
C     ..
C     .. External Subroutines ..
      LOGICAL TEST
      EXTERNAL CINIT,ZCOPY,ZSCAL,ZGEMM,TEST
      INTRINSIC ABS,ZABS
c
C     .. Save statement ..
      SAVE
c ------------------------------------------------------------------------
C     ..



      DO 70 NS = 1,NSHELL
        I = NSH1(NS)
        J = NSH2(NS)
c
        LDIA = 
     +       (DABS(RATOM(1,NS)**2+
     +            RATOM(2,NS)**2+
     +            RATOM(3,NS)**2  ) .LT. 1.0D-6)
c
        DO 10 IU = 1,NSYMAT
c
c --->    GLL = sum(i=1,iumax)(tauvbz * DLL(i) * GS * DLL(i)^T)
c
          IF (IU.EQ.1) THEN
c
c --->      ull(1) is equal to unity matrix
c
            CALL ZCOPY(LMMAXSO*LMMAXSO,GS(1,1,1,NS),1,GLL,1)
            CALL ZSCAL(LMMAXSO*LMMAXSO,TAUVBZ,GLL,1)
C            CALL ZGEMM('N','T',LMAXSQ,LMAXSQ,LMAXSQ,CONE,TPG,LMAXSQ,
C    +           ULL(1,1,IU),LMAXSQ,CZERO,GLL,LMAXSQ)
          ELSE
c
c --->      tpg = tauvbz * DLL * GS
c                         N
          CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,TAUVBZ,
     +                DSYMLL(1,1,IU),LMMAXSO,GS(1,1,IU,NS),
     +                LMMAXSO,CZERO,TPG,LMMAXSO)

c         WRITE(36,*) "IU",IU 
c         DO LM1=1,LMMAXSO
c           DO LM2=1,LMMAXSO
c             WRITE(36,"((2I5),(2e17.9))") LM1,LM2,TPG(LM2,LM1)
c           END DO
c         END DO
           
c
c --->    GLL = GLL + TPG * DLL(i)^T
c                           T
            CALL ZGEMM('N','T',LMMAXSO,LMMAXSO,LMMAXSO,CONE,TPG,LMMAXSO,
     +           DSYMLL(1,1,IU),LMMAXSO,CONE,GLL,LMMAXSO)
          END IF
c     
 10     CONTINUE                    ! IU = 1,NSYMAT
c
c --->  XC = TINVLL(I) * GLL
c
        CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,TINVLL(1,1,I),
     +       LMMAXSO,GLL,LMMAXSO,CZERO,XC,LMMAXSO)
        
c         WRITE(36,*) "XC"
c         WRITE(36,*) "NS",NS
c         DO LM1=1,LMMAXSO
c           DO LM2=1,LMMAXSO
c             WRITE(36,"((2I5),(2e17.9))") LM1,LM2,XC(LM2,LM1)
c           END DO
c         END DO
           
           
        IF (LDIA) THEN
c           WRITE(6,*) "LDIA true"
c     
c --->    GLL = -TINVLL - TINVLL * GLL* TINVLL
c     
          CALL ZCOPY(LMMAXSO**2,TINVLL(1,1,I),1,GLL,1)
          CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,-CONE,XC,LMMAXSO,
     +         TINVLL(1,1,I),LMMAXSO,-CONE,GLL,LMMAXSO)

c          WRITE(37,*) "GLL"
c          WRITE(37,*) "NS",NS
c         DO LM1=1,LMMAXSO
c           DO LM2=1,LMMAXSO
c             WRITE(37,"((2I5),(2e17.9))") LM1,LM2,GLL(LM2,LM1)
c           END DO
c         END DO
           
        ELSE                        ! (LDIA)
c           WRITE(6,*) "LDIA false"
c     
c --->    GLL =  - TINVLL(I) * GLL* TINVLL(J)
c
          CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,-CONE,XC,LMMAXSO,
     +         TINVLL(1,1,J),LMMAXSO,CZERO,GLL,LMMAXSO)
          
        END IF                      ! (LDIA)
c
c --->  GLL0 = GLL/RFCTOR
c
        DO 30 LM1 = 1,LMMAXSO
          DO 20 LM2 = 1,LMMAXSO
            GLL0(LM2,LM1,NS) = GLL(LM2,LM1)/RFCTOR
 20       CONTINUE
 30     CONTINUE

 70   CONTINUE                      ! NS = 1,NSHELL

      RETURN
 9000 format(2f18.10)

      END




