C ************************************************************************
      SUBROUTINE WFINT0_SO(CDER,DDER,QZLM,QZEKDR,PZEKDR,VNSPLL,LSM,
     +                HSOFAC,IRMIN,NSRA,LMMAXD,NSPIN)
C ************************************************************************
c      determines the integrands CDER, DDER or ADER, BDER in the
c        integral equations for the non-spherical wavefunctions from
c        the non-spherical contributions of the potential vinsPLL.
c        (This subroutine is used in zeroth order Born approximation,
c         otherwise subroutine WFINT must be used)
c      R. Zeller      Aug. 1994
C ************************************************************************
C     .. Parameters ..
      include 'inc.p'
      INTEGER LMMAXD,NSPIN
c      PARAMETER (LMMAXD= (LMAXD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
C     ..
C     .. Scalar Arguments ..
      INTEGER IRMIN,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX CDER(2*LMMAXD,2*LMMAXD,IRMIN:IRMD),
     +               DDER(2*LMMAXD,2*LMMAXD,IRMIN:IRMD),
     +               PZEKDR(LMMAXD,IRMIN:IRMD,2,NSPIND),
     +               QZEKDR(LMMAXD,IRMIN:IRMD,2,NSPIND),
     +               QZLM(LMMAXD,IRMIN:IRMD,2,NSPIND),
     +               LSM(2*LMMAXD,2*LMMAXD)
      DOUBLE PRECISION HSOFAC(IRMIN:IRMD)
      DOUBLE PRECISION VNSPLL(2*LMMAXD,2*LMMAXD,IRMIN:IRMD)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX V1
      INTEGER IR,LM1,LM2,LM1MOD,LM2MOD
C     ..
c     OPEN (UNIT=31,file="QZEKDR",form="formatted")
c     OPEN (UNIT=32,file="PZEKDR",form="formatted")
c     OPEN (UNIT=33,file="PZLM",form="formatted")

      DO 50 IR = IRMIN,IRMD

c test Swantje      
c       DO LM2 = 1,LMMAXD
c         WRITE(31,"((2I5),(4e17.9))") IR,LM2,
c    +           QZEKDR(LM2,IR,1)
c         WRITE(32,"((2I5),(4e17.9))") IR,LM2,
c    +           PZEKDR(LM2,IR,1)
c         WRITE(33,"((2I5),(4e17.9))") IR,LM2,
c    +           QZLM(LM2,IR,1)
c       END DO

        DO LM2 = 1,LMMAXD
          DO LM1 = 1,LMMAXD
c            write(44,*) IR,LM1,LM2,VNSPLL(LM1,LM2,IR)
            V1 = (VNSPLL(LM1,LM2,IR)+LSM(LM1,LM2)*HSOFAC(IR))
     +                                    *QZLM(LM2,IR,1,1)
            CDER(LM1,LM2,IR) = QZEKDR(LM1,IR,1,1)*V1
            DDER(LM1,LM2,IR) = PZEKDR(LM1,IR,1,1)*V1
          END DO
        END DO

        DO LM2 = LMMAXD+1,2*LMMAXD
          LM2MOD=MOD(LM2-1,LMMAXD)+1
          DO LM1 = LMMAXD+1,2*LMMAXD
            LM1MOD=MOD(LM1-1,LMMAXD)+1
            V1 = (VNSPLL(LM1,LM2,IR)+LSM(LM1,LM2)*
     +                           HSOFAC(IR))*QZLM(LM2MOD,IR,1,NSPIN)
            CDER(LM1,LM2,IR) = QZEKDR(LM1MOD,IR,1,NSPIN)*V1
            DDER(LM1,LM2,IR) = PZEKDR(LM1MOD,IR,1,NSPIN)*V1
          END DO
        END DO

        DO LM2 = 1,LMMAXD
          DO LM1 = LMMAXD+1,2*LMMAXD
            V1 = (LSM(LM1,LM2)*HSOFAC(IR))*QZLM(LM2,IR,1,1)
            CDER(LM1,LM2,IR) = QZEKDR(LM1-LMMAXD,IR,1,NSPIN)*V1
            DDER(LM1,LM2,IR) = PZEKDR(LM1-LMMAXD,IR,1,NSPIN)*V1
          END DO
        END DO

        DO LM2 = LMMAXD+1,2*LMMAXD
          DO LM1 = 1,LMMAXD
            V1 = (LSM(LM1,LM2)*HSOFAC(IR))*QZLM(LM2-LMMAXD,IR,1,NSPIN)
            CDER(LM1,LM2,IR) = QZEKDR(LM1,IR,1,1)*V1
            DDER(LM1,LM2,IR) = PZEKDR(LM1,IR,1,1)*V1
          END DO
        END DO


        IF (NSRA.EQ.2) THEN

          DO LM2 = 1,LMMAXD
            DO LM1 = 1,LMMAXD
              V1 = (VNSPLL(LM1,LM2,IR)+LSM(LM1,LM2)*HSOFAC(IR))
     +                                    *QZLM(LM2,IR,2,1)
              CDER(LM1,LM2,IR) = CDER(LM1,LM2,IR) +QZEKDR(LM1,IR,2,1)*V1
              DDER(LM1,LM2,IR) = DDER(LM1,LM2,IR) +PZEKDR(LM1,IR,2,1)*V1
            END DO   
          END DO   

          DO LM2 = LMMAXD+1,2*LMMAXD
            LM2MOD=MOD(LM2-1,LMMAXD)+1
            DO LM1 = LMMAXD+1,2*LMMAXD
              LM1MOD=MOD(LM1-1,LMMAXD)+1
              V1 = (VNSPLL(LM1,LM2,IR)+LSM(LM1,LM2)*
     +                           HSOFAC(IR))*QZLM(LM2MOD,IR,2,NSPIN)
              CDER(LM1,LM2,IR) = CDER(LM1,LM2,IR)
     +                              + QZEKDR(LM1MOD,IR,2,NSPIN)*V1
              DDER(LM1,LM2,IR) = DDER(LM1,LM2,IR) + 
     +                                PZEKDR(LM1MOD,IR,2,NSPIN)*V1
            END DO   
          END DO   

          DO LM2 = LMMAXD+1,2*LMMAXD
            DO LM1 = 1,LMMAXD
              V1 = (LSM(LM1,LM2)*HSOFAC(IR))*QZLM(LM2-LMMAXD,IR,1,NSPIN)
              CDER(LM1,LM2,IR) = CDER(LM1,LM2,IR) +QZEKDR(LM1,IR,2,1)*V1
              DDER(LM1,LM2,IR) = DDER(LM1,LM2,IR) +PZEKDR(LM1,IR,2,1)*V1
            END DO   
          END DO   

          DO LM2 = 1,LMMAXD
            DO LM1 = LMMAXD+1,2*LMMAXD
              V1 = (LSM(LM1,LM2)*HSOFAC(IR))*QZLM(LM2,IR,1,1)
              CDER(LM1,LM2,IR) = CDER(LM1,LM2,IR) + 
     +                                 QZEKDR(LM1-LMMAXD,IR,2,NSPIN)*V1
              DDER(LM1,LM2,IR) = DDER(LM1,LM2,IR) + 
     +                                 PZEKDR(LM1-LMMAXD,IR,2,NSPIN)*V1
            END DO   
          END DO   

        END IF
      DO LM2=1,2*LMMAXD
       Do LM1=1,2*LMMAXD
c       write(44,*) IR,LM1,LM2,CDER(LM1,LM2,IR),DDER(LM1,LM2,IR)
       ENDDO
      ENDDO

   50 CONTINUE

c     CLOSE(31)
c     CLOSE(32)
c     CLOSE(33)

      END
