MODULE mod_rhocore_kkrimp

  CONTAINS
  !------------------------------------------------------------------------------------
  !> Summary: Driver for the density from core electrons.
  !> Category: core-electrons, KKRimp
  !------------------------------------------------------------------------------------
  SUBROUTINE RHOCORE(ISPIN,NSPIN,IATOM,NATOM,CELL,VISP,ZATOM,CORESTATE,DENSITY,CONFIG,ITSCF)

    !                          DRDI,R,VISP,A,B,ZATOM,IRCUT, &
    !                          RHOC,ECORE,NCORE,LCORE, &
    !                           )
    ! ,CSCL,
    !                        VTREL,BTREL,RMREL,
    !      &                   DRDIREL,R2DRDIREL,ZREL,JWSREL,IRSHIFT,ECOREREL,
    !      &                   NKCORE,KAPCORE)
    !
    use nrtype
    use mod_rhocoreint
    use mod_config, only: config_testflag
    use type_cell
    use type_config
    use type_corestate
    use type_density
    use mod_types, only: t_inc
    IMPLICIT NONE

    !       INCLUDE 'inc.p'
    ! INTEGER         :: NSRA
    INTEGER,intent(in)              ::  ISPIN
    INTEGER,intent(in)              ::  NSPIN
    INTEGER,intent(in)              ::  IATOM
    INTEGER,intent(in)              ::  NATOM
    ! type(cell_type)                 :: cell
    TYPE(CELL_TYPE)                 ::  CELL
    real(kind=DP),intent(in)        ::  VISP(CELL%NRMAX)
    real(kind=DP),intent(in)        ::  ZATOM
    ! real(kind=DP)                   ::  RHOC(CELL%NRMAX,NSPIN)
    type(corestate_type)            ::  corestate                ! derived type containing all shape information 
    type(density_type)              ::  density                
    type(config_type)              ::  config
    integer                         :: ITSCF


    ! real(kind=DP)   :: ECORE(20))
    ! INTEGER         :: NCORE
    ! INTEGER         :: LCORE(20)

    !     ..!
    !     .. Scalar Arguments ..
    !       DOUBLE PRECISION A,B,ZATOM
    ! !       INTEGER JWSREL,ZREL,IRSHIFT
    !       INTEGER IATOM,ISPIN,NCORE,NSPIN,NSRA
    ! C     ..
    ! C     .. Array Arguments ..
    !       DOUBLE PRECISION DRDI(IRMD),ECORE(20)),&
    !                        R(IRMD),RHOC(IRMD,2),&
    !                        VISP(IRMD)
    !       INTEGER IRCUT(0:IPAND),&
    !               LCORE(20)
    ! ===================================================================
    !  RELATIVISTIC TREATMENT OF CORE ELECTRONS   July/2002
    !  SEE ROUTINE <DRVCORE> FOR A SHORT DESCRIPTION OF THE VARIABLES
    ! C
    !       DOUBLE PRECISION ECOREREL(KREL*20+(1-KREL),2)
    !       INTEGER NKCORE(20),KAPCORE(20*2)
    !       DOUBLE PRECISION CSCL(KREL*LMAXD+1)
    !       DOUBLE PRECISION VTREL(IRMD*KREL+(1-KREL))
    !       DOUBLE PRECISION BTREL(IRMD*KREL+(1-KREL))
    !       DOUBLE PRECISION DRDIREL(IRMD*KREL+(1-KREL)),
    !      &                 R2DRDIREL(IRMD*KREL+(1-KREL)),
    !      &                 RMREL(IRMD*KREL+(1-KREL))
    ! C ===================================================================
    !     ..
    !     .. Local Scalars ..
          DOUBLE PRECISION QC,QC1,RMAX
          INTEGER IPR,NR
          SAVE QC
    !     ..
    !     .. External Subroutines ..
    !       EXTERNAL COREL,DRVCORE
    ! --------------------------------------------------------------
    !     ipr=0 : do not write state dependent information
    !     ipr=1 : write something
    !     ipr=2 : write all (for debugging)
    ! --------------------------------------------------------------
    !          write(*,*) 'sdf'
    if (iatom==1 .and. ispin==1) then
      if (t_inc%i_write>0) write(1337,*) '----------------------------------------------------'
      if (t_inc%i_write>0) write(1337,*) '-------           MODULE RHOCORE             -------'
      if (t_inc%i_write>0) write(1337,*) '----------------------------------------------------'


    end if
    if (ispin==1 .and. ITSCF==1) then
      allocate(density%rhoc(CELL%NRMAX,nspin))
    end if

          IPR = 0
    !     
          IF ( ISPIN.EQ.1 ) QC = 0.0D0

          NR = CELL%NRCUT(1)
    !       NR = IRCUT(1)
          RMAX = CELL%RMESH(NR)
    !       RMAX = R(NR)

    !=======================================================================
    ! non/scalar-relativistic OR relativistic
    !
    !       IF ( KREL.EQ.0 ) THEN
    !     

    ! CELL%LOGPARAMS()
             CALL RHOCOREINT(config%NSRA,IPR,IATOM,DENSITY%RHOC(:,ISPIN),VISP, & 
                             corestate%ECORE(:,ISPIN),corestate%LCORE(:,ISPIN),corestate%NCORE, &
                             cell%DRMESHDI, ZATOM, QC1, & 
                             CELL%LOGPARAMS(1),CELL%LOGPARAMS(2), & 
                             ISPIN,NSPIN,NR,RMAX,CELL%NRMAX)
    !          CALL COREL(NSRA,IPR,IATOM,RHOC(1,ISPIN),VISP,ECORE,LCORE,NCORE,
    !      +        DRDI,ZATOM,QC1,A,B,ISPIN,NSPIN,NR,RMAX,IRMD)


    !          write(*,*) ispin
             IF (IPR.NE.0) WRITE (*,FMT=99001) IATOM
             QC = QC + QC1
             IF (ISPIN.EQ.NSPIN) THEN 
               if (t_inc%i_write>0) write(1337,*) 'ATOM',IATOM
               if (t_inc%i_write>0) WRITE (1337,FMT=99002) ZATOM,QC
             END IF
    !=======================================================================
    !       ELSE
    ! C=======================================================================
    !          CALL DRVCORE(IPR,IATOM,LCORE,NCORE,CSCL,VTREL,BTREL,RMREL,A,B,
    !      &                DRDIREL,R2DRDIREL,ZREL,JWSREL,IRSHIFT,RHOC,
    !      &                ECOREREL,NKCORE,KAPCORE,ECORE,LMAXD,IRMD)
    !       END IF
    ! C
    ! C non/scalar-relativistic OR relativistic
    !=======================================================================
    !       RETURN
    99001 FORMAT (1x,5('*'),' core-relaxation for ',i3,'th cell', &
                  ' was done ',5('*'))
    99002 FORMAT (4X,'nuclear charge  ',F10.6,9X,'core charge =   ',F10.6)

    if (config_testflag('writercore')) then
      if (ispin==1 .and. iatom==1) then
        open(unit=10000,file='test_rcore')
      end if
      write(10000,*) '#iatom',iatom
      write(10000,*) '#ispin',ispin
      write(10000,'(50000g24.16)') density%rhoc(:,ispin)
      if (ispin==nspin .and. iatom==natom) then
        close(10000)
      end if
    end if

  END SUBROUTINE RHOCORE

END MODULE mod_rhocore_kkrimp
