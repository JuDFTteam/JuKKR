MODULE MOD_RITES
  CONTAINS
    !---------------------------------------------------------------------
    !> Summary: Writes the potential file
    !> Category: input-output, KKRimp, potential
    !>   this subroutine stores in 'ifile' the necessary results
    !>   (potentials e.t.c.) to start self-consistency iterations
    !>    modified for the full potential case - if ins .gt. 0 there
    !>    is written a different potential card
    !>    if the sum of absolute values of an lm component of vins (non
    !>    spher. potential) is less than the given rms error qbound this
    !>    component will not be stored .
    !>        (see to subroutine start , where most of the arrays are
    !>      described)
    !>                            modified by b. drittler  aug. 1988
    !--------------------------------------------------------------------
    SUBROUTINE RITES(NSPIN,NATOM,ZATOM,ALAT,NRMAXD,LMAXD,INS, &
                       QBOUND,EFERMI,VBC,CELL,CORESTATE,LMAXATOM,VPOT_OUT)
! C     .. Parameters ..
!       include 'inc.p'
      USE TYPE_CELL
      USE TYPE_CORESTATE
      USE MOD_CONFIG, only : config_testflag
      use mod_version_info
      IMPLICIT NONE

      INTEGER              :: NSPIN
      INTEGER              :: NATOM
      REAL*8               :: ZATOM(NATOM)
      REAL*8               :: ALAT
      INTEGER              :: NRMAXD
      INTEGER              :: LMAXD
!       INTEGER              :: KXC
      INTEGER              :: INS
      REAL*8               :: QBOUND
      REAL*8               :: EFERMI
      REAL*8               :: VBC(2)
      TYPE(CELL_TYPE)      :: CELL(NATOM)
      TYPE(CORESTATE_TYPE) :: CORESTATE(NATOM)
      INTEGER              :: LMAXATOM(NATOM)
      real*8               :: vpot_out(:,:,:,:) !thetas(iri,nfund,*),
!       real*8               :: vpot_out(nrmaxd,(2*lmaxd+1)**2,nspin,natom) !thetas(iri,nfund,*),



!       PARAMETER (LMPOTD= (LPOTD+1)**2)
!       INTEGER IRMIND
!       PARAMETER (IRMIND=IRMD-IRNSD)
! C     ..
! C     .. Scalar Arguments ..
!       DOUBLE PRECISION ALAT,QBOUND
      INTEGER IFILE !,INS,KSHAPE,KXC,LPOT,NATPS,NATYP,NSPIN
! C     ..
! C     .. Array Arguments ..
!       DOUBLE PRECISION A(*),B(*),DRDI(IRMD,*),ECORE(20,*),EFERMI, !(2), 22.5,2000
!      +                 R(IRMD,*),RMT(*),RMTNEW(*),RWS(*),VBC(2),
!      +                 VINS(IRMIND:IRMD,LMPOTD,*),VM2Z(IRMD,*),Z(*)
! C ===================================================================
!       DOUBLE PRECISION ECOREREL(KREL*20+(1-KREL),
!      +                          2*NATYPD)  ! relativistic core energies
!       INTEGER NKCORE(20,NATYPD),KAPCORE(20,2*NATYPD)
! C ===================================================================
!       INTEGER IRC(*),IRNS(*),IRWS(*),ITITLE(20,*),LCORE(20,*),NCORE(*)
      CHARACTER*124 TXC(10)
! C     ..
! C     .. Local Scalars ..
      DOUBLE PRECISION A1,B1,RMAX,RMT1,RMTNW1,RV,SIGN,SUM,Z1
      INTEGER ICORE,INEW,IP,IR,IRMIN,IRNS1,ISAVE,LM,LMNR, &
              LMPOT,NCORE1,NR,IATOM,ISPIN,LPOT
! C     ..
! C     .. Local Arrays ..
      DOUBLE PRECISION DRADI(NRMAXD),ECORE1(20),RA(NRMAXD),VM2ZA(NRMAXD)
      INTEGER LCORE1(20)
! C     ..
! C     .. Intrinsic Functions ..
      INTRINSIC SQRT
! C     ..
! c -------------------------------------------------------------------
      ISAVE = 1
      INEW  = 1
      IFILE = 324249
      TXC(1) = ' Morruzi,Janak,Williams  #serial: ' // serialnr
      TXC(2) = ' von Barth,Hedin         #serial: ' // serialnr
      TXC(3) = ' Vosko,Wilk,Nusair       #serial: ' // serialnr
      TXC(4) = ' GGA PW91                #serial: ' // serialnr

      OPEN(UNIT=IFILE,FILE='out_potential')
  if ( config_testflag('write_density') ) then
      OPEN(UNIT=4357324,FILE='test_potential')
  end if
! c
! c
      DO IATOM = 1,NATOM
        DO ISPIN = 1,NSPIN
          LPOT = 2*LMAXATOM(IATOM)
          LMPOT = (2*LMAXATOM(IATOM)+1)**2

          if ( config_testflag('write_density') ) then
            DO LM=1,(2*LMAXATOM(IATOM)+1)**2
              WRITE(4357324,'(5000g25.16)') VPOT_OUT(:,LM,ISPIN,IATOM)
            END DO
          end if !config_testflag('write_density')
          IF (ISPIN.EQ.NSPIN) THEN
            SIGN = 1.0D0
          ELSE
            SIGN = -1.0D0
          END IF

          IP = NSPIN* (IATOM-1) + ISPIN

!           RMT1 = RMT(IATOM)
          RMT1 = CELL(IATOM)%RMT

!           RMTNW1 = RMTNEW(IATOM)
          RMTNW1 = CELL(IATOM)%RMT !???

          Z1 = ZATOM(IATOM)

!           RMAX = RWS(IATOM)
          RMAX = CELL(IATOM)%RMAX

          IF (INS.EQ.0) THEN
            NR = CELL(IATOM)%NRMAX ! NR = IRWS(IATOM)
          ELSE
            NR = CELL(IATOM)%NRMAX !IRC(IATOM)        !    NR = IRC(IATOM)
          END IF

          IRNS1 = CELL(IATOM)%NRNS !          IRNS1 = IRNS(IATOM)
          IRMIN = NR - IRNS1

          if (ins==1) then
            IF (IRMIN/=CELL(IATOM)%NRMIN_NS) stop '[rites] error IRMIN/=NRMIN_NS'
          end if

          A1 = CELL(IATOM)%LOGPARAMS(1) !A(IATOM)
          B1 = CELL(IATOM)%LOGPARAMS(2) !B(IATOM)


!           CORESTATE needs to be changed from (IATOM) to (ISPIN,IATOM)
          NCORE1 = CORESTATE(IATOM)%NCORE !(IP)
! c
          DO IR = 1,NR
            RA(IR) = CELL(IATOM)%RMESH(IR) !R(IR,IATOM)
            DRADI(IR) = CELL(IATOM)%DRMESHDI(IR) !DRDI(IR,IATOM)
! c
! c--->       store only lm=1 component of the potential
! c
            VM2ZA(IR) = VPOT_OUT(IR,1,ISPIN,IATOM) !VM2Z(IR,IP)
          END DO !J
! C
          WRITE (IFILE,FMT=9000) CELL(IATOM)%VPOT_NAME(ISPIN),TXC(CELL(IATOM)%KXC+1)
          WRITE (IFILE,*) RMT1,ALAT,RMTNW1
          !WRITE (IFILE,FMT=9010) RMT1,ALAT,RMTNW1
          WRITE (IFILE,FMT=9020) Z1,RMAX,EFERMI,VBC(ISPIN)
          WRITE (IFILE,FMT=9030) NR,A1,B1,NCORE1,INEW
! C
          IF (NCORE1.GE.1) THEN
! C
!             IF (KREL.EQ.0) THEN 
               DO ICORE = 1,NCORE1
                 LCORE1(ICORE) = CORESTATE(IATOM)%LCORE(ICORE,ISPIN)  !LCORE(ICORE,IP)
                 ECORE1(ICORE) = CORESTATE(IATOM)%ECORE(ICORE,ISPIN) !ECORE(ICORE,IP)
               END DO
               WRITE (IFILE,FMT=9040) (LCORE1(ICORE),&
                    ECORE1(ICORE),ICORE=1,NCORE1)
!             ELSE
!               DO J = 1,NCORE1
!                  LCORE1(J) = LCORE(J,IP)
!                  ECORE2(J,1) = ECOREREL(J,2*IH-1)
!                  ECORE2(J,2) = ECOREREL(J,2*IH)
!               END DO
! C
! C --> independent of spin, the \mu-averaged relativistic core energies
! C     are written out for \kappa = -l-1,l 
! C     format compatible with the non-(scalar) relativistic mode
! C     however, the next read in has no meaning for the REL core-solver
! C     a detailed output of the core energies is supplied by < CORE >
! C
!               DO ICORE=1,NCORE1
!                  WRITE (IFILE,FMT=9041) LCORE1(ICORE),(
!      +                ECORE2(ICORE,I+1),
!      +                TXTL(LCORE1(ICORE)),
!      +                TXTK(IABS(KAPCORE(ICORE,2*IH-1+I))),
!      +                I=0,NKCORE(ICORE,IH)-1)
!               END DO
!             END IF
! C
          END IF
! c

          IF (INS.EQ.0) THEN
! c
! c--->       store only the spherically averaged potential 
! c           (in mt or as - case) 
! c           this is done always for the host
! c
            IF (INEW.EQ.0) THEN
              WRITE (IFILE,FMT=9050) &
                   (RA(IR),DRADI(IR),VM2ZA(IR),IR=1,NR)
            ELSE
              WRITE (IFILE,FMT=9051) (VM2ZA(IR),IR=1,NR)
            END IF

          ELSE
! c
! c--->     store the full potential , but the non spherical contribution
! c         only from irns1 up to irws1 ;
! c         remember that the lm = 1 contribution is multiplied
! c         by a factor 1/sqrt(4 pi)
! c
            WRITE (IFILE,FMT=9060) NR,IRNS1,LMPOT,ISAVE
            WRITE (IFILE,FMT=9070) (VM2ZA(IR),IR=1,NR)
            IF (LPOT.GT.0) THEN
              LMNR = 1
              DO LM = 2,LMPOT
                SUM = 0.0D0
                DO IR = IRMIN,NR
                  RV = VPOT_OUT(IR,LM,ISPIN,IATOM)*RA(IR)
                  SUM = SUM + RV*RV*DRADI(IR)
                END DO

                IF (SQRT(SUM).GT.QBOUND) THEN
                  LMNR = LMNR + 1
                  WRITE (IFILE,FMT=9060) LM
                  WRITE (IFILE,FMT=9070) (VPOT_OUT(IR,LM,ISPIN,IATOM),IR=IRMIN,NR)
                END IF

              END DO
! c
! c--->         write a one to mark the end
! c
              IF (LMNR.LT.LMPOT) WRITE (IFILE,FMT=9060) ISAVE
            END IF

          END IF

        END DO !ISPIN
      END DO !IATOM

      CLOSE(IFILE)

 9000 FORMAT (a28,6x,'  exc:',a124,3x,a10)
 9010 FORMAT (3f12.8)
 9020 FORMAT (f10.5,/,f10.5,2f15.10)
 9030 FORMAT (i3,/,2d15.8,/,2i2)
 9040 FORMAT (i5,1p,d20.11)
 9041 FORMAT (i5,2(1p,d20.11,2x,A1,A3))
 9050 FORMAT (1p,2d15.6,1p,d15.8)
 9051 FORMAT (1p,4d20.12)
 9060 FORMAT (10i5)
 9070 FORMAT (1p,4d20.13)
      END SUBROUTINE
      END MODULE MOD_RITES
