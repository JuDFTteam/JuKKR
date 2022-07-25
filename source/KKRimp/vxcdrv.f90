!------------------------------------------------------------------------------------
!> Summary: Driver for the exchange-correlation potential and energy calculation
!> Author: 
!> Driver for the exchange-correlation potential and energy calculation. It wraps
!> all the different exchange-correlation potentials and make sure to call the 
!> appropriate subroutines depending on the type of exchange correlation potential
!> indicated in the `inputcard`
!------------------------------------------------------------------------------------
      MODULE MOD_VXCDRV
      CONTAINS
!   call vxcdrv(energyparts%exc,config%kte,config%kxc,nspin,natom,density, & 
!               vpot_out, cell,config%kshape,gauntshape, shapefun,lmaxd, & 
!               (2*lmaxd), (2*lmaxd+1)**2, (4*lmaxd+1)**2, cell(1)%nrmaxd, lmaxatom)
  !-------------------------------------------------------------------------------
  !> Summary: Driver for the exchange-correlation potential and energy calculation
  !> Author:
  !> Category: xc-potential, KKRimp
  !> Deprecated: False 
  !> Driver for the exchange-correlation potential and energy calculation. It wraps
  !> all the different exchange-correlation potentials and make sure to call the 
  !> appropriate subroutines depending on the type of exchange correlation potential
  !> indicated in the `inputcard`
  !-------------------------------------------------------------------------------
      SUBROUTINE VXCDRV(EXC,KTE,NSPIN,NATOM,DENSITY,VONS, &
                        CELL,KSHAPE,gauntshape,&
                        SHAPEFUN,LMAXD,LPOTD,LMPOTD,LMXSPD,nrmaxd,lmaxatom,ins)
      USE TYPE_CELL
      USE TYPE_DENSITY
      USE TYPE_SHAPEFUN
      USE TYPE_GAUNTSHAPE
      USE MOD_EXCHANGECORRELATION
      use mod_vxcspo, only: vxcspo
      use mod_vxclm, only: vxclm
      use mod_vxcgga, only: vxcgga
      use global_variables, only: ipand, irmd, ngshd, irid, nfund
      IMPLICIT NONE

!       INCLUDE 'inc.p'
!interface
       DOUBLE PRECISION       :: EXC(0:LPOTD,NATOM)
       INTEGER                :: KTE
!        INTEGER                :: KXC
       INTEGER                :: LPOT
       INTEGER                :: NSPIN
       INTEGER                :: NATOM
       type(density_type)     :: density(natom)
       DOUBLE PRECISION       :: VONS(nrmaxd,LMPOTD,NSPIN,NATOM)
       TYPE(CELL_TYPE)        ::  CELL(NATOM)
       INTEGER                :: kshape
       TYPE(GAUNTSHAPE_TYPE)  :: GAUNTSHAPE(LMAXD)
       type(shapefun_type)    :: shapefun(natom)
       INTEGER                ::  LMAXD
       INTEGER                :: LPOTD
       INTEGER                :: LMPOTD
       INTEGER                :: LMXSPD
       INTEGER                :: nrmaxd
       INTEGER                :: lmaxatom(natom)
       INTEGER                :: INS
!local
       INTEGER                :: LMAX
       INTEGER                :: LMPOT
       INTEGER, parameter     :: IJD=434 !LMXSPD

      DOUBLE PRECISION        :: DYLMF1(IJD,LMPOTD),DYLMF2(IJD,LMPOTD), &
                                 DYLMT1(IJD,LMPOTD),DYLMT2(IJD,LMPOTD), &
                                 DYLMTF(IJD,LMPOTD),RHO2IAT(nrmaxd,LMPOTD,2), &
                                 RIJ(IJD,3), RIJ_NOGGA(IJD,3),&
                                 THET(IJD),& 
                                 WTYR(IJD,LMPOTD),WTYR_NOGGA(IJD,LMPOTD), &
                                 YLM(IJD,LMPOTD),&
                                 YR(IJD,LMPOTD),YR_NOGGA(IJD,LMPOTD)
      INTEGER                 :: IFUNMIAT(LMXSPD),LMSPIAT(LMXSPD)
      INTEGER                 :: IATOM,ISPIN,LMX1,LM

!       IF (KXC.LT.3) THEN
        CALL SPHERE_NOGGA(LPOTD,YR_NOGGA,WTYR_NOGGA,RIJ_NOGGA,IJD)
!         write(144,*) LPOTD,YR_NOGGA,WTYR_NOGGA,RIJ_NOGGA,IJD
!       ELSE
!         write(154,*)LPOTD,YR,WTYR,RIJ,IJD,LMPOTD,THET,YLM,DYLMT1, &
!                         DYLMT2,DYLMF1,DYLMF2,DYLMTF
        CALL SPHERE_GGA(LPOTD,YR,WTYR,RIJ,IJD,LMPOTD,THET,YLM,DYLMT1, &
                        DYLMT2,DYLMF1,DYLMF2,DYLMTF)
!         write(155,*) LPOTD,YR,WTYR,RIJ,IJD,LMPOTD,THET,YLM,DYLMT1, &
!                         DYLMT2,DYLMF1,DYLMF2,DYLMTF
!       END IF
!         write(145,*) LPOTD,YR_NOGGA,WTYR_NOGGA,RIJ_NOGGA,IJD

      DO IATOM = 1,NATOM
        LMAX =  LMAXATOM(IATOM)
        LPOT =2*LMAXATOM(IATOM)
        LMPOT=(LPOT+1)**2

        IFUNMIAT=0
        LMSPIAT=0

        if (ins/=0) then
          DO LMX1 = 1,(4*LMAXATOM(IATOM)+1)**2
            IFUNMIAT(LMX1) = SHAPEFUN(IATOM)%LM2INDEX(LMX1)
            LMSPIAT(LMX1) =  SHAPEFUN(IATOM)%LMUSED(LMX1)
          END DO
        END IF

        RHO2IAT=0.0D0
        DO ISPIN=1, NSPIN
          DO LM=1,(2*LMAXATOM(IATOM)+1)**2
            RHO2IAT(:CELL(IATOM)%NRMAX,LM,ISPIN)=DENSITY(IATOM)%RHO2NS(:,LM,ISPIN)
          END DO !LM
        END DO !ISPIN
! *************  old **********************************
!         CALL DCOPY(nrmaxd*LMPOTD,DENSITY(IATOM)%RHO2NS(1,1,1),1,RHO2IAT(1,1,1),1)
!         IF (NSPIN.EQ.2) THEN
!           CALL DCOPY(nrmaxd*LMPOTD,DENSITY(IATOM)%RHO2NS(1,1,IATOM,2),1,RHO2IAT(1,1,2),1)
!         END IF

!           write(244,*) EXC,KTE,KXC,LPOT,NSPIN,IATOM
!           write(245,*) RHO2IAT
!                      VONS(:,:,:,IATOM),CELL(IATOM)%RMESH,CELL(IATOM)%DRMESHDI,&
!                      CELL(IATOM)%NRMAX,CELL(IATOM)%NRCUT,CELL(IATOM)%NPAN,&
!                      KSHAPE,GAUNTSHAPE(LMAX)%GSH,GAUNTSHAPE(LMAX)%ILM,GAUNTSHAPE(LMAX)%IMAXSH,IFUNMIAT,SHAPEFUN(IATOM)%THETAS,&
!                      YR_NOGGA,WTYR_NOGGA,IJD,LMSPIAT,&
!                      LMPOTD,2*LMAXD,LMXSPD,nrmaxd,SHAPEFUN(IATOM)%NLMSHAPED,SHAPEFUN(IATOM)%NRSHAPED,GAUNTSHAPE(LMAX)%NGSHD,&
!                      CELL(IATOM)%NPAND



        ! these values will be used in vxclm from global_variables module, so we set them here for each atom
        irmd = nrmaxd
        nfund = SHAPEFUN(IATOM)%NLMSHAPED
        irid = SHAPEFUN(IATOM)%NRSHAPED
        ngshd = GAUNTSHAPE(LMAX)%NGSHD
        ipand = CELL(IATOM)%NPAND

        IF (CELL(IATOM)%KXC.LT.3) THEN
          

          CALL VXCLM(EXC,KTE,CELL(IATOM)%KXC,LPOT,NSPIN,IATOM,RHO2IAT, &
                     VONS(:,:,:,IATOM),CELL(IATOM)%RMESH,CELL(IATOM)%DRMESHDI,&
                     CELL(IATOM)%NRMAX,CELL(IATOM)%NRCUT,CELL(IATOM)%NPAN,&
                     KSHAPE,GAUNTSHAPE(LMAX)%GSH,GAUNTSHAPE(LMAX)%ILM,GAUNTSHAPE(LMAX)%IMAXSH,IFUNMIAT,SHAPEFUN(IATOM)%THETAS,&
                     YR_NOGGA,WTYR_NOGGA,IJD,LMSPIAT)!,&
!                    lmpotd, lpotd, lmxspd, irmd , nfund                   , irid                   , ngshd, 
                     !LMPOTD,2*LMAXD,LMXSPD,nrmaxd,SHAPEFUN(IATOM)%NLMSHAPED,SHAPEFUN(IATOM)%NRSHAPED,GAUNTSHAPE(LMAX)%NGSHD,&
!                    ipand
                     !CELL(IATOM)%NPAND) 


!      +                 LMPOTD,LPOTD,LMXSPD,IRMD,NFUND,IRID,NGSHD,IPAND)
        ELSE
! c
! c GGA EX-COR POTENTIAL
! c
          CALL VXCGGA(EXC,KTE,CELL(IATOM)%KXC,LPOT,NSPIN,IATOM,RHO2IAT, &
                      VONS(:,:,:,IATOM),CELL(IATOM)%RMESH,CELL(IATOM)%DRMESHDI,CELL(IATOM)%LOGPARAMS(1), &
                      CELL(IATOM)%NRMAX,CELL(IATOM)%NRCUT,CELL(IATOM)%NPAN, &
                      KSHAPE,GAUNTSHAPE(LMAX)%GSH,GAUNTSHAPE(LMAX)%ILM,GAUNTSHAPE(LMAX)%IMAXSH,IFUNMIAT,SHAPEFUN(IATOM)%THETAS, &
                      !YR,WTYR,IJD,LMSPIAT,THET,YLM,DYLMT1,DYLMT2, &
                      WTYR,IJD,LMSPIAT,THET,YLM,DYLMT1,DYLMT2, &
                      DYLMF1,DYLMF2,DYLMTF)!, &
                      !LMPOTD,2*LMAXD,LMXSPD,nrmaxd,SHAPEFUN(IATOM)%NLMSHAPED,SHAPEFUN(IATOM)%NRSHAPED,GAUNTSHAPE(LMAX)%NGSHD,&
                      !CELL(IATOM)%NPAND) 
!      +                 LMPOTD,LPOTD,LMXSPD,IRMD,NFUND,IRID,NGSHD,IPAND)
        END IF
      END DO



      END SUBROUTINE VXCDRV
      END MODULE MOD_VXCDRV
