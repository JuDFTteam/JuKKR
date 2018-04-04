      SUBROUTINE CALC_RHO_LL_SS(LMMAX,RLL,IRCUT,IPAN,ICELL,THETAS,
     +                CLEB,ICLEB,IEND,IFUNM,LMSP,IRWS,DRDI,DENS)
c     +                LM1,LM2)

      IMPLICIT NONE

      include 'inc.p'
      INTEGER          LMMAXD
      PARAMETER        (LMMAXD= (LMAXD+1)**2)
      INTEGER          LMPOTD
      PARAMETER        (LMPOTD= (LPOTD+1)**2)
      INTEGER          IRMIND
      PARAMETER        (IRMIND=IRMD-IRNSD)
C     ..
C     .. Scalar Arguments ..
      INTEGER          IEND,LMMAX,IRWS!,LM1,LM2
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX   RLL(IRMD,LMMAXD,LMMAXD),   ! non-sph. eigen states of single pot 
     +                 DENS
      DOUBLE PRECISION CLEB(*),
     +                 THETAS(IRID,NFUND,*),
     +                 DRDI(IRMD)                            ! derivative dr/di
      INTEGER          ICLEB(NCLEB,4),IFUNM(NATYPD,LMPOTD),
     +                 LMSP(NATYPD,*),IRCUT(0:IPAND),IPAN,
     +                 ICELL,IFUN


c local variables

      DOUBLE PRECISION             ::   C0LL
      DOUBLE COMPLEX               ::   CLT 
      DOUBLE COMPLEX, ALLOCATABLE  ::   RSP(:),RGES(:)!,DENS(:,:,:)
c     +                                  RGES_W(:,:,:,:),
c     +                                  DENS_GESAMT(:,:),
c     +                                  DENS_GESAMT_I1(:,:,:)
      integer                      ::   LM1P,LM2P,LM3P,IR,J,I
      INTEGER                      ::   IRCUTM(0:IPAND)

C     ..
c  ---> first calculate only the spherically symmetric contribution
c       (for all points r; if r>r_MT (or IR> IRMIN),the density has to
c       multiplied with the shape functions...
c
c---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
c
      ALLOCATE(RGES(IRMD))
      ALLOCATE(RSP(IRMD))

c      WRITE(6,*) "In rho ll"

      C0LL = 1.0d0/DSQRT(16.0D0*DATAN(1.0D0))
      RSP=0d0
      RGES=0d0
      
      DO LM1P = 1,LMMAX
        DO IR=1,IRMD
          RSP(IR)=RSP(IR)+RLL(IR,LM1P,LM1P)
        END DO
      END DO
      
      DO 60 IR = 1,IRCUT(IPAN)
        RGES(IR) = RSP(IR)
   60 CONTINUE
           
      IF (IPAN.GT.1) THEN
        DO 100 IR = IRCUT(1)+1,IRCUT(IPAN)
c          WRITE(56,"((I5),(2e17.9))") IR,THETAS(IR-IRCUT(1),1,ICELL)
          RGES(IR) = RSP(IR)*C0LL*THETAS(IR-IRCUT(1),1,ICELL)
  100   CONTINUE
      END IF

c      STOP " "
!      WRITE(6,*) "IRCUT(1)",IRCUT(1)
!      WRITE(6,*) "IRCUT(IPAN)",IRCUT(IPAN)
!      WRITE(6,*) "IRCUT(IPAN)-IRMIND",IRCUT(IPAN)-IRMIND
!      WRITE(6,*) "IRMIND",IRMIND

      DO 160 J = 1,IEND
        LM1P = ICLEB(J,1)
        LM2P = ICLEB(J,2)
        LM3P = ICLEB(J,3)
        CLT = CLEB(J)
 
c---> calculate the non spherically symmetric contribution
 
        IF (IPAN.GT.1 .AND. LMSP(ICELL,LM3P).GT.0) THEN
          IFUN = IFUNM(ICELL,LM3P)
c          WRITE(156,*) "IFUN",IFUN
          IF (LM1P == LM2P ) THEN
            DO 150 IR = IRCUT(1)+1,IRCUT(IPAN)
             RGES(IR) = RGES(IR)+RLL(IR,LM2P,LM1P)*
     +             CLEB(J)*THETAS(IR-IRCUT(1),IFUN,ICELL)
  150       CONTINUE
          ELSE
            DO IR = IRCUT(1)+1,IRCUT(IPAN)
              RGES(IR) = RGES(IR)+
     +             CLEB(J)*THETAS(IR-IRCUT(1),IFUN,ICELL)*
     +          (RLL(IR,LM2P,LM1P)+RLL(IR,LM1P,LM2P))
            END DO
          END IF
        END IF

  160 CONTINUE

      IF (IPAN.EQ.1) THEN
        IRCUTM(0) = 0
        IRCUTM(1) = IRWS
      ELSE
        DO 30 I = 0,IPAN
          IRCUTM(I) = IRCUT(I)
   30   CONTINUE
      END IF

      CALL CSIMPK(RGES(:),DENS,IPAN,IRCUTM,DRDI)

      DEALLOCATE(RGES)
      DEALLOCATE(RSP)

      END SUBROUTINE
