      SUBROUTINE WFINT_SO(QNS,CDER,DDER,QZEKDR,PZEKDR,VNSPLL,LSM,HSOFAC,
     +                       NSRA,IRMIN,IRMD,LMMAXD,NSPIN)
c     Implicit None
c-----------------------------------------------------------------------
c      determines the integrands CDER, DDER or ADER, BDER in the
c        integral equations for the non-spherical wavefunctions from
c        the non-spherical contributions of the potential vinsPLL.
c
c      R. Zeller      Aug. 1994
c
c     changed by Swantje, 21.06.10 for SOC:
c     splitting of the arguments in real and imaginary parts and use
c     of the DGEMM routine is no longer correct since the spin-orbit
c     Hamiltonian is complex.
c-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER IRMD,IRMIN,LMMAXD,NSRA,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX CDER(2*LMMAXD,2*LMMAXD,IRMIN:IRMD),
     +               DDER(2*LMMAXD,2*LMMAXD,IRMIN:IRMD),
     +               PZEKDR(LMMAXD,IRMIN:IRMD,2,NSPIN),
     +               QNS(2*LMMAXD,2*LMMAXD,IRMIN:IRMD,2),
     +               QZEKDR(LMMAXD,IRMIN:IRMD,2,NSPIN),
     +               LSM(2*LMMAXD,2*LMMAXD),
     +               HSO(2*LMMAXD,2*LMMAXD)
      DOUBLE PRECISION HSOFAC(IRMIN:IRMD)
      DOUBLE PRECISION VNSPLL(2*LMMAXD,2*LMMAXD,IRMIN:IRMD)
C     ..
C     .. Local Scalars ..
      INTEGER IR,LM1,LM2,LM1MOD,LM2MOD
C     ..
C     .. External Subroutines ..
      EXTERNAL DGEMM,ZGEMM
C     ..
C     .. Local Arrays ..
c     DOUBLE PRECISION VTQNSI(2*LMMAXD,2*LMMAXD),
c    +                 VTQNSR(2*LMMAXD,2*LMMAXD)
      DOUBLE COMPLEX   VTQNSI(2*LMMAXD,2*LMMAXD),
     +                 VTQNSR(2*LMMAXD,2*LMMAXD),
     +                   QNSI(2*LMMAXD,2*LMMAXD),
     +                   QNSR(2*LMMAXD,2*LMMAXD)
      DOUBLE COMPLEX    VTQNS(2*LMMAXD,2*LMMAXD)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DCMPLX,DIMAG,DBLE
C     ..
c
c      WRITE(6,*) "IRMIN",IRMIN
      DO 90 IR = IRMIN,IRMD
        VTQNSR=0d0
        VTQNSI=0d0
        QNSR=0d0
        QNSI=0d0
        DO 20 LM2 = 1,2*LMMAXD
          DO 10 LM1 = 1,2*LMMAXD
c            IF (HSOFAC(IR) .NE. 0 ) STOP "HSOFAC"
            HSO(LM1,LM2)=HSOFAC(IR)*LSM(LM1,LM2)+
     +                                    VNSPLL(LM1,LM2,IR)
            QNSR(LM1,LM2) = (QNS(LM1,LM2,IR,1))
c            QNSR(LM1,LM2) = REAL(QNS(LM1,LM2,IR,1))
c            QNSI(LM1,LM2) = DIMAG(QNS(LM1,LM2,IR,1))
   10     CONTINUE
   20   CONTINUE

        CALL ZGEMM('N','N',2*LMMAXD,2*LMMAXD,2*LMMAXD,1.D0,
     +          HSO(1:2*LMMAXD,1:2*LMMAXD),2*LMMAXD,
     +             QNSR(1:2*LMMAXD,1:2*LMMAXD),
     +          2*LMMAXD,0.D0,VTQNSR(1:2*LMMAXD,1:2*LMMAXD),2*LMMAXD)

c       CALL ZGEMM('N','N',2*LMMAXD,2*LMMAXD,2*LMMAXD,1.D0,
c    +          HSO(1:2*LMMAXD,1:2*LMMAXD),2*LMMAXD,
c    +             QNSR(1:2*LMMAXD,1:2*LMMAXD),
c    +          2*LMMAXD,0.D0,VTQNSR(1:2*LMMAXD,1:2*LMMAXD),2*LMMAXD)

c       CALL ZGEMM('N','N',2*LMMAXD,2*LMMAXD,2*LMMAXD,1.D0,
c    +          HSO(1:2*LMMAXD,1:2*LMMAXD),2*LMMAXD,
c    +             QNSI(1:2*LMMAXD,1:2*LMMAXD),
c    +          2*LMMAXD,0.D0,VTQNSI(1:2*LMMAXD,1:2*LMMAXD),2*LMMAXD)

        DO LM1 = 1,2*LMMAXD
          LM1MOD=MOD(LM1-1,LMMAXD)+1
          DO LM2 = 1,2*LMMAXD
            LM2MOD=MOD(LM2-1,LMMAXD)+1
            CDER(LM1,LM2,IR) = QZEKDR(LM1MOD,IR,1,1)*
     +         VTQNSR(LM1,LM2)
            DDER(LM1,LM2,IR) = PZEKDR(LM1MOD,IR,1,1)*
     +         VTQNSR(LM1,LM2)
c           CDER(LM1,LM2,IR) = QZEKDR(LM1MOD,IR,1,1)*
c    +         DCMPLX(REAL(VTQNSR(LM1,LM2)),REAL(VTQNSI(LM1,LM2)))
c           DDER(LM1,LM2,IR) = PZEKDR(LM1MOD,IR,1,1)*
c    +         DCMPLX(REAL(VTQNSR(LM1,LM2)),REAL(VTQNSI(LM1,LM2)))
          END DO
        END DO

        IF (NSRA.EQ.2) THEN
          VTQNSR=0d0
          VTQNSI=0d0
          QNSR=0d0
          QNSI=0d0
          DO LM2 = 1,2*LMMAXD
            DO LM1 = 1,2*LMMAXD
c              IF (HSOFAC(IR) .NE. 0 ) STOP "HSOFAC"
              HSO(LM1,LM2)=HSOFAC(IR)*LSM(LM1,LM2)+
     +                                      VNSPLL(LM1,LM2,IR)
              QNSR(LM1,LM2) = QNS(LM1,LM2,IR,2)
c              QNSR(LM1,LM2) = REAL(QNS(LM1,LM2,IR,2))
c              QNSI(LM1,LM2) = DIMAG(QNS(LM1,LM2,IR,2))
            END DO    
          END DO    

          CALL ZGEMM('N','N',2*LMMAXD,2*LMMAXD,2*LMMAXD,1.D0,
     +            HSO(1:2*LMMAXD,1:2*LMMAXD),2*LMMAXD,
     +               QNSR(1:2*LMMAXD,1:2*LMMAXD),
     +            2*LMMAXD,0.D0,VTQNSR(1:2*LMMAXD,1:2*LMMAXD),2*LMMAXD)
         

c         CALL ZGEMM('N','N',2*LMMAXD,2*LMMAXD,2*LMMAXD,1.D0,
c    +            HSO(1:2*LMMAXD,1:2*LMMAXD),2*LMMAXD,
c    +               QNSR(1:2*LMMAXD,1:2*LMMAXD),
c    +            2*LMMAXD,0.D0,VTQNSR(1:2*LMMAXD,1:2*LMMAXD),2*LMMAXD)
c        
c         CALL ZGEMM('N','N',2*LMMAXD,2*LMMAXD,2*LMMAXD,1.D0,
c    +            HSO(1:2*LMMAXD,1:2*LMMAXD),2*LMMAXD,
c    +               QNSI(1:2*LMMAXD,1:2*LMMAXD),
c    +            2*LMMAXD,0.D0,VTQNSI(1:2*LMMAXD,1:2*LMMAXD),2*LMMAXD)


          DO LM1 = 1,2*LMMAXD
            LM1MOD=MOD(LM1-1,LMMAXD)+1
            DO LM2 = 1,2*LMMAXD
              LM2MOD=MOD(LM2-1,LMMAXD)+1
              CDER(LM1,LM2,IR) = CDER(LM1,LM2,IR)+
     +           QZEKDR(LM1MOD,IR,2,1)*
     +           VTQNSR(LM1,LM2)
              DDER(LM1,LM2,IR) = DDER(LM1,LM2,IR)+ 
     +           PZEKDR(LM1MOD,IR,2,1)*
     +           VTQNSR(LM1,LM2)
c             CDER(LM1,LM2,IR) = CDER(LM1,LM2,IR)+
c    +           QZEKDR(LM1MOD,IR,2,1)*
c    +           DCMPLX(REAL(VTQNSR(LM1,LM2)),REAL(VTQNSI(LM1,LM2)))
c             DDER(LM1,LM2,IR) = DDER(LM1,LM2,IR)+ 
c    +           PZEKDR(LM1MOD,IR,2,1)*
c    +           DCMPLX(REAL(VTQNSR(LM1,LM2)),REAL(VTQNSI(LM1,LM2)))
            END DO
          END DO
c         DO 80 LM2 = 1,2*LMMAXD
c           DO 70 LM1 = 1,LMMAXD
c             CDER(LM1,LM2,IR) = CDER(LM1,LM2,IR) +
c    +                           QZEKDR(LM1,IR,2,1)*VTQNS(LM1,LM2)
c             DDER(LM1,LM2,IR) = DDER(LM1,LM2,IR) +
c    +                           PZEKDR(LM1,IR,2,1)*VTQNS(LM1,LM2)
c  70       CONTINUE
c           DO LM1 = LMMAXD+1,2*LMMAXD
c             LM1MOD=MOD(LM1-1,LMMAXD)+1
c             CDER(LM1,LM2,IR) = CDER(LM1,LM2,IR) +
c    +                      QZEKDR(LM1MOD,IR,2,NSPIN)*
c    +         DCMPLX(REAL(VTQNSR(LM1,LM2)),REAL(VTQNSI(LM1,LM2)))
c             DDER(LM1,LM2,IR) = DDER(LM1,LM2,IR) +
c    +                      PZEKDR(LM1MOD,IR,2,NSPIN)*
c    +         DCMPLX(REAL(VTQNSR(LM1,LM2)),REAL(VTQNSI(LM1,LM2)))
c           END DO
c  80     CONTINUE

        END IF

c       IF (NSRA.EQ.2) THEN

c         CALL ZGEMM('N','N',2*LMMAXD,2*LMMAXD,2*LMMAXD,1.D0,HSO,
c    +             2*LMMAXD,QNS(:,:,IR,2),2*LMMAXD,0.D0,VTQNS,2*LMMAXD)

c         DO 80 LM2 = 1,2*LMMAXD
c           DO 70 LM1 = 1,LMMAXD
c             CDER(LM1,LM2,IR) = CDER(LM1,LM2,IR) +
c    +                           QZEKDR(LM1,IR,2,1)*VTQNS(LM1,LM2)
c             DDER(LM1,LM2,IR) = DDER(LM1,LM2,IR) +
c    +                           PZEKDR(LM1,IR,2,1)*VTQNS(LM1,LM2)
c  70       CONTINUE
c           DO LM1 = LMMAXD+1,2*LMMAXD
c             LM1MOD=MOD(LM1-1,LMMAXD)+1
c             CDER(LM1,LM2,IR) = CDER(LM1,LM2,IR) +
c    +                      QZEKDR(LM1MOD,IR,2,NSPIN)*VTQNS(LM1,LM2)
c             DDER(LM1,LM2,IR) = DDER(LM1,LM2,IR) +
c    +                      PZEKDR(LM1MOD,IR,2,NSPIN)*VTQNS(LM1,LM2)
c           END DO
c  80     CONTINUE

c       END IF

   90 CONTINUE

      END
