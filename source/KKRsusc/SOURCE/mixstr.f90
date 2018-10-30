      MODULE MOD_MIXSTR
      CONTAINS
! c 13.10.95 ***************************************************************
      SUBROUTINE MIXSTR(RMSAVQ,RMSAVM,INS,NATOM,LMAXATOM,LMAXD, &
                       NSPIN,ITC, &
                       MIXING,FCM,VPOT, VPOT_OUT,CELL,NRMAXD,lmpotin &
                       )
   USE TYPE_CELL, only: cell_type
   USE MOD_CONFIG, only: CONFIG_RUNFLAG
   IMPLICIT NONE
! c ************************************************************************
! C     .. Parameters ..
!interface
      DOUBLE PRECISION             ::  RMSAVM,RMSAVQ
      INTEGER                      ::  INS
      INTEGER                      ::  NATOM
      INTEGER                      ::  LMAXATOM(NATOM)
      INTEGER                      ::  LMAXD
      INTEGER                      ::  NSPIN
      INTEGER                      ::  ITC
      real*8                       ::  MIXING
      real*8                       ::  FCM ! ????
      real*8                       ::  vpot    (nrmaxd,       lmpotin,nspin,natom) !thetas(iri,nfund,*),
      real*8                       ::  vpot_out(nrmaxd,(2*lmaxd+1)**2,nspin,natom) !thetas(iri,nfund,*),
      TYPE(CELL_TYPE)              ::  CELL(NATOM)
      INTEGER                      ::  NRMAXD
      INTEGER                      :: lmpotin
!local
      real*8                       ::  PI,FPI,RFPI
      INTEGER                      ::  IATOM, LVAL
      INTEGER                      ::  IPF ! some write out file
      INTEGER                      ::  LMPOT
      INTEGER                      ::  LPOT
      INTEGER                      ::  IH,IHP1
      real*8                       ::  FAC,RMSERM,RMSERQ,VMN,VNM,VNP,VOLDM,VOLDP,VPN
      INTEGER                      ::  IRC1,IRMIN1,IR,LM
      INTEGER,SAVE                 ::  IATOMMAXDISPLAY=10
! C     ..
! C     .. Intrinsic Functions ..
      INTRINSIC MOD,REAL,SQRT
      IPF=6
      PI = 4.D0*ATAN(1.D0)
      FPI = 4.D0*PI
      RFPI = SQRT(FPI)

! C     ..
      RMSAVQ = 0.0D0
      RMSAVM = 0.0D0

      IF (CONFIG_RUNFLAG('rmsdisplayall')) IATOMMAXDISPLAY=1000000

! c
! c---> final construction of the potentials
! c     attention : the spherical averaged potential is the lm=1
! c                     component of vons times sqrt(4 pi).
! c
! c     first mixing scheme : straight mixing
! c---> determination of the root mean sqare error
! c
!       NATOM = 0.0D0
      DO IATOM = 1,NATOM
        LVAL=LMAXATOM(IATOM)
        LMPOT=(2*LMAXATOM(IATOM)+1)**2

!         I = IATOM - NATREF
!         NATOM = NATOM + DBLE(NSHELL(I))*CONC(I)

        IF (NSPIN.EQ.2) THEN
          IH = 1
          IHP1 = 2
        ELSE
          IH = 1
          IHP1 = 1
        END IF

        IRC1 = CELL(IATOM)%NRMAX ! IRC(IATOM)
!         write(*,*) IRC1
        RMSERQ = 0.0D0
        RMSERM = 0.0D0
        FAC = 0.5D0/RFPI
! c
        DO IR = 1,IRC1
          VNP = FAC* (vpot_out(IR,1,IH,IATOM)+vpot_out(IR,1,IHP1,IATOM))
          VNM = FAC* (vpot_out(IR,1,IH,IATOM)-vpot_out(IR,1,IHP1,IATOM))
          VOLDP = 0.5D0* (vpot(IR,1,IH,IATOM)+vpot(IR,1,IHP1,IATOM))
          VOLDM = 0.5D0* (vpot(IR,1,IH,IATOM)-vpot(IR,1,IHP1,IATOM))
          RMSERQ = RMSERQ + 2.0D0*REAL(1+MOD(IR,2))*CELL(IATOM)%RMESH(IR)*CELL(IATOM)%RMESH(IR)* &
                   CELL(IATOM)%DRMESHDI(IR)* (VNP-VOLDP)* (VNP-VOLDP)

          RMSERM = RMSERM + 2.0D0*REAL(1+MOD(IR,2))*CELL(IATOM)%RMESH(IR)*CELL(IATOM)%RMESH(IR)* &
                   CELL(IATOM)%DRMESHDI(IR)* (VNM-VOLDM)* (VNM-VOLDM)
          VPN = VOLDP + MIXING* (VNP-VOLDP)
          VMN = VOLDM + FCM*MIXING* (VNM-VOLDM)
          IF (NSPIN==2) vpot_out(IR,1,2,IATOM) = VPN - VMN
          vpot_out(IR,1,1,IATOM) = VPN + VMN
        END DO
!        write(*,*) CELL(IATOM)%RMESH(IRC1)**3
        RMSERQ = RMSERQ/ (CELL(IATOM)%RMESH(IRC1)**3)
        RMSERM = RMSERM/ (CELL(IATOM)%RMESH(IRC1)**3)
        RMSAVQ = RMSAVQ + RMSERQ !*NSHELL(I)*CONC(I)
        RMSAVM = RMSAVM + RMSERM !*NSHELL(I)*CONC(I)

        IF (NSPIN.EQ.2) THEN
          WRITE (1337,FMT=9000) IATOM,SQRT(RMSERQ),SQRT(RMSERM)
          IF (NATOM<=IATOMMAXDISPLAY) THEN
            WRITE (   *,FMT=9000) IATOM,SQRT(RMSERQ),SQRT(RMSERM)
          END IF
        ELSE
          WRITE (1337,FMT=9020) IATOM,SQRT(RMSERQ)
          IF (NATOM<=IATOMMAXDISPLAY) THEN
            WRITE (   *,FMT=9020) IATOM,SQRT(RMSERQ)
          END IF
        END IF

        IF (INS.NE.0) THEN
          RMSERQ = 0.0D0
          RMSERM = 0.0D0
          IRMIN1 = CELL(IATOM)%NRMIN_NS !IRMIN(IATOM)
          DO LM = 2,LMPOT
            DO IR = IRMIN1,IRC1
              VNP = 0.5D0* (vpot_out(IR,LM,IH,IATOM)+vpot_out(IR,LM,IHP1,IATOM))
              VNM = 0.5D0* (vpot_out(IR,LM,IH,IATOM)-vpot_out(IR,LM,IHP1,IATOM))
              VOLDP = 0.5D0* (vpot(IR,LM,IH,IATOM)+vpot(IR,LM,IHP1,IATOM))
              VOLDM = 0.5D0* (vpot(IR,LM,IH,IATOM)-vpot(IR,LM,IHP1,IATOM))
              RMSERQ = RMSERQ + 2.0D0*REAL(1+MOD(IR,2))*CELL(IATOM)%RMESH(IR)*CELL(IATOM)%RMESH(IR)* &
                       CELL(IATOM)%DRMESHDI(IR)* (VNP-VOLDP)* (VNP-VOLDP)
              RMSERM = RMSERM + 2.0D0*REAL(1+MOD(IR,2))*CELL(IATOM)%RMESH(IR)*CELL(IATOM)%RMESH(IR)* &
                       CELL(IATOM)%DRMESHDI(IR)* (VNM-VOLDM)* (VNM-VOLDM)
              VPN = VOLDP + MIXING* (VNP-VOLDP)
              VMN = VOLDM + FCM*MIXING* (VNM-VOLDM)
              vpot_out(IR,LM,IHP1,IATOM) = VPN - VMN
              vpot_out(IR,LM,IH,IATOM) = VPN + VMN
            END DO
          END DO
          RMSERQ = RMSERQ/ (CELL(IATOM)%RMESH(IRC1)**3)/FPI
          RMSERM = RMSERM/ (CELL(IATOM)%RMESH(IRC1)**3)/FPI
          RMSAVQ = RMSAVQ + RMSERQ !*NSHELL(I)*CONC(I)
          RMSAVM = RMSAVM + RMSERM !*NSHELL(I)*CONC(I)

          IF (NSPIN.EQ.2) THEN
            WRITE (1337,FMT=9010) IATOM,SQRT(RMSERQ),SQRT(RMSERM)
            IF (NATOM<=IATOMMAXDISPLAY) THEN
              WRITE (   *,FMT=9010) IATOM,SQRT(RMSERQ),SQRT(RMSERM)
            END IF
          ELSE
            WRITE (1337,FMT=9030) IATOM,SQRT(RMSERQ)
            IF (NATOM<=IATOMMAXDISPLAY) THEN
              WRITE (   *,FMT=9030) IATOM,SQRT(RMSERQ)
            END IF
          END IF

        END IF

      END DO ! IATOM

      RMSAVQ = SQRT(RMSAVQ/NATOM)
      RMSAVM = SQRT(RMSAVM/NATOM)

      WRITE(6,'(79(1H-),/)')
      IF (NSPIN.EQ.2) THEN
        WRITE (1337,FMT=9040) ITC,RMSAVQ,RMSAVM
        WRITE (   *,FMT=9040) ITC,RMSAVQ,RMSAVM
      ELSE
        WRITE (1337,FMT=9050) ITC,RMSAVQ
        WRITE (   *,FMT=9050) ITC,RMSAVQ
      END IF
!       WRITE(1337,'(79(1H-))')
!       WRITE(   *,'(79(1H-))')

 9000 FORMAT (5x,' rms-error for atom',i3,1x,':','v+ + v- = ',1p,d11.4, &
            2x,',',2x,'v+ - v- = ',1p,d11.4)
 9010 FORMAT (5x,' rms-error non spherical contribution for atom ',i3, &
            1x,':','v+ + v- = ',1p,d11.4,02x,',',2x,'v+ - v- = ',1p, &
            d11.4)
 9020 FORMAT (5x,' rms-error for atom',i3,1x,':','v+ + v- = ',1p,d11.4)
 9030 FORMAT (5x,' rms-error non spherical contribution for atom ',i3, &
            1x,':','v+ + v- = ',1p,d11.4)
 9040 FORMAT ('      ITERATION',I4,' average rms-error : v+ + v- = ', &
            1p,d11.4,/,39x,' v+ - v- = ',1p,d11.4)
 9050 FORMAT ('      ITERATION',I4,' average rms-error : v+ + v- = ', &
            1p,d11.4)
      END SUBROUTINE
      END MODULE 