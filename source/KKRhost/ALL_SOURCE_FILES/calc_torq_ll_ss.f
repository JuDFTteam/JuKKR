c     This subroutine computes a matrix that is the basis for constructing
c     the KKR representation of the torque operator. It is adapted from the 
c     CALC_RHO_LL_SS subroutine, but the spin dependent part, i.e., the exhange
c     field, replaces the shape function in the integration.
c
c                                     Guillaume Geranton, September 2014
      SUBROUTINE CALC_TORQ_LL_SS(LMMAX,RLL,IRCUT,IPAN,ICELL,
     +                CLEB,ICLEB,IEND,IFUNM,LMSP,IRWS,DRDI,DENS,
     +                VISP,NSPIN,IATOM,VINS,IRMIN)
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
      INTEGER          IEND,LMMAX,IRWS,NSPIN,IATOM,IRMIN!,LM1,LM2
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX   RLL(IRMD,LMMAXD,LMMAXD),   ! non-sph. eigen states of single pot 
     +                 DENS
      DOUBLE PRECISION CLEB(*),
     +                 DRDI(IRMD),                            ! derivative dr/di
     +                 VISP(IRMD,*), !              spherical part of the potential
     +                 VINS(IRMIND:IRMD,LMPOTD,*) ! non-sph. part of the potential
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

      EXTERNAL TEST
      LOGICAL TEST

C     ..
c  ---> first calculate only the spherically symmetric contribution
c       (for all points r; if r>r_MT (or IR> IRMIN),the density has to
c       multiplied with the shape functions...
c
c---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
c
      ALLOCATE(RGES(IRMD))
      ALLOCATE(RSP(IRMD))

      C0LL = 1.0d0/DSQRT(16.0D0*DATAN(1.0D0))
      RSP=0d0
      RGES=0d0

c     Compute spherical contribution to the torque (LM=1)
c     Sph. potential has to be multiplied by sqrt(4 PI) !
      DO LM1P = 1,LMMAX
       DO IR=1,IRMD
         RSP(IR)=RSP(IR)+RLL(IR,LM1P,LM1P)*C0LL*(-1)*
     &   (VISP(IR,NSPIN*(IATOM-1)+2)-VISP(IR,NSPIN*(IATOM-1)+1))*0.5*
     &   DSQRT(16.0D0*DATAN(1.0D0)) 

       END DO
      END DO

      DO 60 IR = 1,IRMD
c cut contributions from outside the MT if recquired
        IF ((TEST('ONLYMT  ') .eq. .TRUE.) .and. (IR > IRCUT(1))) THEN
          RGES(IR) = 0
        ELSE
          RGES(IR) = RSP(IR)
        END IF
   60 CONTINUE

      IF (TEST('ONLYSPH ') .eq. .FALSE.) THEN
        DO 160 J = 1,IEND
          LM1P = ICLEB(J,1)
          LM2P = ICLEB(J,2)
          LM3P = ICLEB(J,3)   ! always >= 2 here
          CLT = CLEB(J)
 
c--->   calculate the non spherically symmetric contribution 
         IF (IPAN.GT.1 .AND. LMSP(ICELL,LM3P).GT.0) THEN
           IFUN = IFUNM(ICELL,LM3P)
           IF (LM1P == LM2P ) THEN
c             DO 150 IR = IRCUT(1)+1,IRCUT(IPAN)
             DO 150 IR = IRMIN,IRCUT(IPAN)
              RGES(IR) = RGES(IR)+RLL(IR,LM2P,LM1P)*CLEB(J)*(-1)*
     +                   (VINS(IR,LM3P,NSPIN* (IATOM-1) + 2) -
     +                   VINS(IR,LM3P,NSPIN* (IATOM-1) + 1))*0.5
150          CONTINUE
           ELSE
c             DO IR = IRCUT(1)+1,IRCUT(IPAN)
             DO IR = IRMIN,IRCUT(IPAN)
               RGES(IR) = RGES(IR)+
     +                CLEB(J)*(-1)*(VINS(IR,LM3P,NSPIN* (IATOM-1) + 2) -
     +                VINS(IR,LM3P,NSPIN* (IATOM-1) + 1))*0.5*
     +                (RLL(IR,LM2P,LM1P)+RLL(IR,LM1P,LM2P))
             END DO
           END IF
         END IF

  160   CONTINUE
      END IF

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
