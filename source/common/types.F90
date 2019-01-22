!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Module defining necessary types for the MPI communication
!> Author: 
!> Module defining necessary types for the MPI communication, as well as memory management
!> and initialization of the needed arrays.
!------------------------------------------------------------------------------------
module mod_types

  use :: mod_datatypes, only: dp
  use :: mod_constants, only: czero

  implicit none

  public :: t_inc, t_tgmat, t_mpi_c_grid, t_lloyd, t_dtmatjij, t_cpa, t_imp


  !-------------------------------------------------------------------------------
  !> Summary: Type holding single site t-matrix, reference GF and structural GF (gmat) for distribution between 1a, 1b and 1c parts of the code
  !> Author: 
  !> Category: single-site, structural-greensfunction, reference-system, KKRhost 
  !> Deprecated: False 
  !> Type holding single site t-matrix, reference GF and structural GF (gmat)
  !> for distribution between 1a, 1b and 1c parts of the code
  !-------------------------------------------------------------------------------
  type :: type_tgmatices

    ! logical switches to control if matrices are stored in memory or written to files
    logical :: tmat_to_file = .false.
    logical :: gmat_to_file = .false.
    logical :: gref_to_file = .false.

    integer :: nelements = 4       ! 3 arrays in this type, for mpi bcast

    ! allocatable arrays for tmat, gmat and gref
    complex (kind=dp), dimension(:,:,:), allocatable :: tmat !! Single-site t-matrix ! dimensions=LMMAXD, LMMAXD, IREC; IREC= IE+IELAST*(ISPIN-1)+IELAST*NSPIN*(I1-1) ;IE=1,...,IELAST, ISPIN=1,...,NSPIN, I1=1,...,NATYP)
    complex (kind=dp), dimension(:,:,:), allocatable :: gmat !! Structural Greens function ! dimensions=LMMAXD, LMMAXD, IREC; IREC= IQDOS+NQDOS*(IE-1)+NQDOS*IELAST*(ISPIN-1)+IELAST*NSPIN*(I1-1) ;IE=1,...,IELAST, ISPIN=1,...,NSPIN, I1=1,...,NATYP)
    complex (kind=dp), dimension(:,:,:,:), allocatable :: gref !! Reference Greens function ! GINP(NACLSD*LMGF0D,LMGF0D,NCLSD) IREC=IE=1,...,IELAST

  end type type_tgmatices

  !-------------------------------------------------------------------------------
  !> Summary: Type holding CPA information
  !> Author: 
  !> Category: coherent-potential-approx, KKRhost
  !> Deprecated: False 
  !> Type holding CPA information
  !-------------------------------------------------------------------------------
  type :: type_cpa

    ! logical switches to control if matrices are stored in memory or written to files
    logical :: dmatproj_to_file = .false.

    integer :: nelements = 3       ! 2 array in this type, for mpi bcast

    ! allocatable arrays for tmat, gmat and gref
    complex (kind=dp), dimension(:,:,:,:), allocatable :: dmatts ! dimensions=LMMAXD, LMMAXD, NATYP, IREC; IREC= IE+IELAST*(ISPIN-1)+; IE=1,...,IELAST, ISPIN=1,...,NSPIN)
    complex (kind=dp), dimension(:,:,:,:), allocatable :: dtilts ! dimensions=LMMAXD, LMMAXD, NATYP, IREC; IREC= IE+IELAST*(ISPIN-1)+; IE=1,...,IELAST, ISPIN=1,...,NSPIN)
  end type type_cpa


  !-------------------------------------------------------------------------------
  !> Summary: Data type for the derivatives of the t-matrix with respect to changing the non-collinear angles in directions {x,y,z}
  !> Author: 
  !> Category: single-site, KKRhost
  !> Deprecated: False 
  !> Data type for the derivatives of the t-matrix with respect to changing 
  !> the non-collinear angles in directions {x,y,z}
  !-------------------------------------------------------------------------------
  !> data type for the derivatives of the t-matrix with respect to changing the non-collinear angles in directions {x,y,z}
  type :: type_dtmatjijdij

    integer :: nelements = 3
    logical :: calculate = .false.
    complex (kind=dp), dimension(:,:,:,:), allocatable :: dtmat_xyz !! Derivatives of the t-matrix with respect to non-collinear angles ! dimensions= LMMAXD, LMMAXD, 3, IELAST;  3={x,y,z}

  end type type_dtmatjijdij


  !-------------------------------------------------------------------------------
  !> Summary: Type holding some array dimensions needed independently of `t_params`
  !> Author: 
  !> Category: initialization, KKRhost
  !> Deprecated: False 
  !> Type holding some array dimensions needed independently of `t_params`
  !-------------------------------------------------------------------------------
  type :: type_inc

    integer :: nparams = 23        ! number of parameters in type_inc, excluding allocatable array KMESH
    integer :: lmmaxd = -1
    integer :: nspin = -1
    integer :: ielast = -1
    integer :: nqdos = -1
    integer :: natyp = -1
    integer :: lmgf0d = -1
    integer :: nclsd = -1
    integer :: naclsmax = -1
    integer :: i_iteration = -1
    integer :: n_iteration = -1
    integer :: mit_bry = 1
    integer :: nshell0 = -1
    integer :: nkmesh = -1
    logical :: newsosol = .false.  !! use new solver for SOC
    logical :: nosoc = .false.     !! use new solver without SOC (test option 'NOSOC   ')
    logical :: deci_out = .false.  !! use deci_out case
    integer :: i_write = 0         !! switch to control if things are written out or not (verbosity levels 0,1,2)
    integer :: i_time = 1          !! switch to control if timing files are written (verbosity levels 0,1,2)
    ! parameters needed for wavefunctions
    integer :: nsra         = -1
    integer :: lmmaxso      = -1
    integer :: irmdnew      = -1
    integer :: kvrel        = -1

    integer, dimension(:), allocatable :: kmesh, kmesh_ie

  end type type_inc

  !-------------------------------------------------------------------------------
  !> Summary: Type holding information on the MPI parallelization scheme
  !> Author: 
  !> Category: initialization, communication, KKRhost
  !> Deprecated: False 
  !> Type holding information on the MPI parallelization scheme
  !-------------------------------------------------------------------------------
  type :: type_mpi_cartesian_grid_info

    integer :: nparams        = 12
    integer :: mympi_comm_ie  = -1
    integer :: mympi_comm_at  = -1
    integer :: myrank_ie      = -1
    integer :: myrank_at      = -1
    integer :: myrank_atcomm  = -1
    integer :: nranks_ie      = -1
    integer :: nranks_at      = -1
    integer :: nranks_atcomm  = -1
    integer :: ntot1          = -1
    integer :: ntot2          = -1
    integer, dimension(2) :: dims = [ -1, -1 ]
    integer, dimension(:), allocatable :: ntot_pt1
    integer, dimension(:), allocatable :: ioff_pt1
    integer, dimension(:), allocatable :: ntot_pt2
    integer, dimension(:), allocatable :: ioff_pt2

  end type type_mpi_cartesian_grid_info


  !-------------------------------------------------------------------------------
  !> Summary: Type holding information needed for lloyd such as derivatives of single site t-matrix, reference GF or the trace of alpha matrix and
  !> Author: 
  !> Category: single-site, reference-system, KKRhost
  !> Deprecated: False 
  !> Type holding information needed for lloyd such as derivatives of single
  !> site t-matrix, reference GF or the trace of alpha matrix and
  !-------------------------------------------------------------------------------
  type :: type_lloyd

    ! logical switches to control if matrices are stored in memory or written to files
    logical :: dtmat_to_file          = .false. ! unit 691
    logical :: tralpha_to_file        = .false. ! unit 692
    logical :: cdos_diff_lly_to_file  = .false. ! unit 701
    logical :: dgref_to_file          = .false. ! unit 681
    logical :: g0tr_to_file           = .false. ! unit 682

    integer :: n1 = 6              ! 5 logicals and 5 arrays this type, for mpi bcast

    ! allocatable arrays
    complex (kind=dp), dimension(:), allocatable :: g0tr ! complex (kind=dp) LLY_G0TR_IE, irec=ie
    complex (kind=dp), dimension(:), allocatable :: tralpha ! complex (kind=dp) TRALPHA, IREC = ie_num + ie_end*(ISPIN-1) + ie_end*NSPIN* (I1-1)
    complex (kind=dp), dimension(:,:), allocatable :: cdos ! complex (kind=dp) CDOS_LLY(IELAST,NSPIND), irec=IE, aalready in dim 1 of cdos!
    complex (kind=dp), dimension(:,:,:), allocatable :: dtmat ! complex (kind=dp) TMAT0(LMMAXD,LMMAXD), IREC = ie_num + ie_end*(ISPIN-1) + ie_end*NSPIN* (I1-1)
    complex (kind=dp), dimension(:,:,:,:), allocatable :: dgref ! complex (kind=dp); ALLOCATE ( DGINP(NACLSMAX*LMGF0D,LMGF0D,NCLS) ), IREC=IE

  end type type_lloyd

  !-------------------------------------------------------------------------------
  !> Summary: Type holding information for impurity potential, needed in GREENIMP mode
  !> Author: 
  !> Category: potential, shape-functions, geometry, KKRhost
  !> Deprecated: False 
  !> Type holding information for impurity potential, needed in GREENIMP mode
  !-------------------------------------------------------------------------------
  type :: type_imp

    integer :: n1       = 12       ! number of scalars for mpi bcast + 2 (for N1,N2)
    integer :: n2       = 16       ! number of arrays for mpi bcast + 1 (for N2)
    integer :: natomimp = -1       ! number of atoms in impurity cluster
    integer :: ihost    = -1       ! number of different host atoms (layer indices)
    !--------------------------------------------------------------------------------
    ! Array dimensions. can be read from t_params 
    !--------------------------------------------------------------------------------
    integer :: irmd   !! Maximum number of radial points
    integer :: irid   !! Shape functions parameters in non-spherical part
    integer :: ipand  !! Number of panels in non-spherical part
    integer :: nfund  !! Shape functions parameters in non-spherical part
    integer :: nspin  !! Counter for spin directions
    integer :: natypd !! Number of kinds of atoms in unit cell
    integer :: irmind !! irmd - irnsd
    integer :: lmpotd !! (lpot+1)**2
    !--------------------------------------------------------------------------------
    ! Array dimensions. can be read from t_params 
    !--------------------------------------------------------------------------------

    ! allocatable arrays
    integer, dimension(:), allocatable :: ipanimp !! Radial mesh, Panel mesh for impurities ! IPANIMP(NATOMIMP)
    integer, dimension(:), allocatable :: irwsimp !! Radial mesh, IRWS for imps             ! IRWSIMP(NATOMIMP)
    integer, dimension(:), allocatable :: hostimp !! Layer index of host atoms              ! HOSTIMP(NATYPD)
    integer, dimension(:), allocatable :: atomimp !! Layer index of imp atoms               ! ATOMIMP(NATOMIMP)
    integer, dimension(:), allocatable :: irminimp !! Radial mesh, IRMIN for imps           ! IRMINIMP(NATOMIMP))
    integer, dimension(:,:), allocatable :: ircutimp !! Radial mesh, RCUT for imps           ! IRCUTIMP(0:IPAND,NATOMIMP)
    real (kind=dp), dimension(:), allocatable :: zimp     !! atom charge of imps,             ! ZIMP(NATOMIMP)
    real (kind=dp), dimension(:), allocatable :: phiimp   !! phi of nonco_angle of impurity   ! PHIIMP(NATOMIMP)
    real (kind=dp), dimension(:), allocatable :: thetaimp !! theta of nonco_angle of impurity ! THETAIMP(NATOMIMP)
    real (kind=dp), dimension(:,:), allocatable :: rimp   !! Rmesh of imps,                   ! RIMP(IRMD,NATOMIMP)
    real (kind=dp), dimension(:,:), allocatable :: rclsimp !! impurity positions(scoef file)  ! RCLSIMP(3,NATOMIMPD)
    real (kind=dp), dimension(:,:), allocatable :: vispimp !! impurity potential              ! VISPIMP(IRMD,NATOMIMP*NSPIN)
    real (kind=dp), dimension(:,:,:), allocatable :: vinsimp !! impurity potential ! VINSIMP(IRMIND:IRMD,LMPOTD,NATOMIMP*NSPIN)
    real (kind=dp), dimension(:,:,:), allocatable :: thetasimp !! shape functions of imps ! THETASIMP(IRID,NFUND,NATOMIMP)
    complex (kind=dp), dimension(:,:,:,:), allocatable :: rllimp !! impurity wavefunctions ! RLL(NVEC*LMMAXSO,LMMAXSO,IRMDNEW(I1))

  end type type_imp

  !-------------------------------------------------------------------------------
  !> Summary: Type holding information for the madelung potentials
  !> Author: Philipp Ruessmann
  !> Category: KKRhost, initialization, geometry
  !> Deprecated: False 
  !> Needed if avmad and abvmad files are not written but kept in memory
  !> 
  !> @warning
  !> There is no MPI communication of this type yet since it is created in main0
  !> part and only used in main2 part which are all done by the master rank. This
  !> needs to be changed if the parallelization is improved int he future.
  !> @endwarning
  !-------------------------------------------------------------------------------
  type :: type_madel

    integer :: n1       = 12       ! number of scalars for mpi bcast + 2 (for N1,N2)
    !--------------------------------------------------------------------------------
    ! Array dimensions. can be read from t_params 
    !--------------------------------------------------------------------------------
    integer :: irmd   !! Maximum number of radial points

    ! allocatable arrays
    real(kind=dp), dimension(:,:,:), allocatable :: avmad !!Structure-dependent matrix, dimension: irec, lmpot x lmpot 
    real(kind=dp), dimension(:,:), allocatable :: bvmad !!Structure-dependent vector, dimension: irec, lmpot

  end type type_madel

  ! save types
  type (type_inc), save :: t_inc
  type (type_tgmatices), save :: t_tgmat
  type (type_mpi_cartesian_grid_info), save :: t_mpi_c_grid
  type (type_lloyd), save :: t_lloyd
  type (type_dtmatjijdij), allocatable, save :: t_dtmatjij(:) ! dimensions I1=1,...,NATYP
  type (type_cpa), save :: t_cpa
  type (type_imp), save :: t_imp
  type (type_madel), save :: t_madel

contains

  !-------------------------------------------------------------------------------
  !> Summary: Subroutine to allocate and initialize arrays of `t_tgmat`
  !> Author: 
  !> Category: initialization, memory-management, structural-greensfunction, reference-system, KKRhost
  !> Deprecated: False 
  !> Subroutine to allocate and initialize arrays of `t_tgmat`
  !-------------------------------------------------------------------------------
  subroutine init_tgmat(t_inc, t_tgmat, t_mpi_c_grid)

    implicit none

    type (type_inc), intent (in) :: t_inc
    type (type_mpi_cartesian_grid_info), intent (in) :: t_mpi_c_grid
    type (type_tgmatices), intent (inout) :: t_tgmat

    integer :: i_stat, nspin

    nspin = t_inc%nspin
    if (t_inc%newsosol .and. .not.t_inc%nosoc) nspin = 1

    if (.not. allocated(t_tgmat%tmat)) then
      if (.not. t_tgmat%tmat_to_file) then
        ! if(nranks.eq.1) then
        ! !allocate tmat(lmmax,lmmax,irec_max) for irec_max=ielast*nspin*natyp
        ! allocate(t_tgmat%tmat(t_inc%LMMAXD,t_inc%LMMAXD,t_inc%IELAST*nspin*t_inc%NATYP), STAT=ierr)
        ! if(ierr/=0) stop 'Problem allocating t_tgmat%tmat'
        ! else
        ! !allocate tmat(lmmax,lmmax,irec_max) for irec_max=iemax_local*nspin*natyp
        ! allocate(t_tgmat%tmat(t_inc%LMMAXD,t_inc%LMMAXD,t_mpi_c_grid%ntot2*nspin*t_inc%NATYP), STAT=ierr)
        ! if(ierr/=0) stop 'Problem allocating t_tgmat%tmat for mpi'
        ! end if
        ! allocate tmat(lmmax,lmmax,irec_max) for irec_max=iemax_local*nspin*natyp
        allocate(t_tgmat%tmat(t_inc%lmmaxd,t_inc%lmmaxd,t_mpi_c_grid%ntot2*nspin*t_inc%natyp),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_tgmat%tmat'
      else
        allocate(t_tgmat%tmat(1,1,1),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_tgmat%tmat'
      end if

      t_tgmat%tmat(:, :, :) = czero
    end if

    if (.not. allocated(t_tgmat%gmat)) then
      if (.not. t_tgmat%gmat_to_file) then
        ! if(nranks.eq.1) then
        ! !allocate gmat(lmmax,lmmax,irec_max) for irec_max=NQDOS*IELAST*NSPIN*NATYP
        ! allocate(t_tgmat%gmat(t_inc%LMMAXD,t_inc%LMMAXD,t_inc%NQDOS*t_inc%IELAST*nspin*t_inc%NSHELL0), STAT=ierr)
        ! if(ierr/=0) stop 'Problem allocating t_tgmat%gmat'
        ! else
        ! !allocate gmat(lmmax,lmmax,irec_max) for irec_max=NQDOS*IEMAX_local*NSPIN*NATYP
        ! allocate(t_tgmat%gmat(t_inc%LMMAXD,t_inc%LMMAXD,t_inc%NQDOS*t_mpi_c_grid%ntot2*nspin*t_inc%NSHELL0), STAT=ierr)
        ! if(ierr/=0) stop 'Problem allocating t_tgmat%gmat for mpi'
        ! end if
        allocate(t_tgmat%gmat(t_inc%lmmaxd,t_inc%lmmaxd,t_inc%nqdos*t_mpi_c_grid%ntot2*nspin*t_inc%nshell0),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_tgmat%gmat'
      else
        allocate(t_tgmat%gmat(1,1,1),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_tgmat%gmat'
      end if
      t_tgmat%gmat(:, :, :) = czero
    end if

    if (.not. allocated(t_tgmat%gref)) then
      if (.not. t_tgmat%gref_to_file) then
        ! if(nranks.eq.1) then
        ! !allocate gref(NACLSMAX*LMGF0D,LMGF0D,NCLSD,irec_max) for irec_max=IELAST
        ! allocate(t_tgmat%gref(t_inc%NACLSMAX*t_inc%LMGF0D,t_inc%LMGF0D,t_inc%NCLSD,t_inc%IELAST), STAT=ierr)
        ! if(ierr/=0) stop 'Problem allocating t_tgmat%gref'
        ! else
        ! !allocate gref(NACLSMAX*LMGF0D,LMGF0D,NCLSD,irec_max) for irec_max=IEMAX_local (=ntot2)
        ! allocate(t_tgmat%gref(t_inc%NACLSMAX*t_inc%LMGF0D,t_inc%LMGF0D,t_inc%NCLSD,t_mpi_c_grid%ntot2), STAT=ierr)
        ! if(ierr/=0) stop 'Problem allocating t_tgmat%gref for mpi'
        ! end if
        allocate(t_tgmat%gref(t_inc%naclsmax*t_inc%lmgf0d,t_inc%lmgf0d,t_inc%nclsd,t_mpi_c_grid%ntot2),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_tgmat%gref'
      else
        allocate(t_tgmat%gref(1,1,1,1),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_tgmat%gref'
      end if
      t_tgmat%gref(:, :, :, :) = czero
    end if

  end subroutine init_tgmat

  !-------------------------------------------------------------------------------
  !> Summary: Subroutine to allocate and initialize arrays of `t_cpa`
  !> Author: 
  !> Category: initialization, memory-management, coherent-potential-approx, KKRhost
  !> Deprecated: False 
  !> Subroutine to allocate and initialize arrays of `t_cpa`
  !-------------------------------------------------------------------------------
  subroutine init_t_cpa(t_inc, t_cpa, nenergy)

    implicit none

    type (type_inc), intent (in) :: t_inc
    integer, intent (in) :: nenergy
    type (type_cpa), intent (inout) :: t_cpa

    integer :: i_stat, nspin

    nspin = t_inc%nspin
    if (t_inc%newsosol .and. .not.t_inc%nosoc) nspin = 1

    if (.not. allocated(t_cpa%dmatts)) then
      if (.not. t_cpa%dmatproj_to_file) then
        ! allocate tmat(lmmax,lmmax,NATYP,irec_max) for irec_max=nenergy*nspin
        allocate(t_cpa%dmatts(t_inc%lmmaxd,t_inc%lmmaxd,t_inc%natyp,nenergy*nspin),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_cpa%dmatts'
      else
        allocate(t_cpa%dmatts(1,1,1,1),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_cpa%dmatts'
      end if
      t_cpa%dmatts(:, :, :, :) = czero
    end if

    if (.not. allocated(t_cpa%dtilts)) then
      if (.not. t_cpa%dmatproj_to_file) then
        ! allocate tmat(lmmax,lmmax,NATYP,irec_max) for irec_max=nenergy*nspin
        allocate(t_cpa%dtilts(t_inc%lmmaxd,t_inc%lmmaxd,t_inc%natyp,nenergy*nspin),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_cpa%dtilts'
      else
        allocate(t_cpa%dtilts(1,1,1,1),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_cpa%dtilts'
      end if
      t_cpa%dtilts(:, :, :, :) = czero
    end if

  end subroutine init_t_cpa

  !-------------------------------------------------------------------------------
  !> Summary: Subroutine to allocate and initialize arrays of `t_dtmatJij`
  !> Author: 
  !> Category: initialization, memory-management, KKRhost
  !> Deprecated: False 
  !> Subroutine to allocate and initialize arrays of `t_dtmatJij`
  !-------------------------------------------------------------------------------
  subroutine init_t_dtmatjij(t_inc, t_dtmatjij)

    implicit none

    type (type_inc), intent (in) :: t_inc
    type (type_dtmatjijdij), intent (inout), allocatable :: t_dtmatjij(:)

    integer :: i_stat

    if (.not. t_inc%newsosol) stop 'in init_t_dtmatJij: should only be called with NEWSOSOL'

    if (.not. allocated(t_dtmatjij)) then
      allocate(t_dtmatjij(t_inc%natyp),stat=i_stat)
      if (i_stat/=0) stop 'Problem allocating t_dtmatjij'
    end if

  end subroutine init_t_dtmatjij

  !-------------------------------------------------------------------------------
  !> Summary: Subroutine to allocate and initialize arrays of `t_dtmatJij_at`
  !> Author: 
  !> Category: initialization, memory-management, KKRhost
  !> Deprecated: False 
  !> Subroutine to allocate and initialize arrays of `t_dtmatJij_at`
  !-------------------------------------------------------------------------------
  subroutine init_t_dtmatjij_at(t_inc, t_mpi_c_grid, t_dtmatjij_at)

    implicit none

    type (type_inc), intent (in) :: t_inc
    type (type_mpi_cartesian_grid_info), intent (in) :: t_mpi_c_grid
    type (type_dtmatjijdij), intent (inout) :: t_dtmatjij_at

    integer :: i_stat

    if (.not. t_inc%newsosol) stop 'in init_t_dtmatJij_single: should only be called with NEWSOSOL'

    if (.not. allocated(t_dtmatjij_at%dtmat_xyz) .and. t_dtmatjij_at%calculate) then
      ! if (.not. t_dtmatJij%dtmat_to_file) then
      ! if(nranks.eq.1) then
      ! !allocate dtmat_xyz(lmmax,lmmax,3,irec) for irec_max=ielast
      ! allocate(t_dtmatJij_at%dtmat_xyz(t_inc%LMMAXD,t_inc%LMMAXD,3,t_inc%IELAST), STAT=ierr)
      ! if(ierr/=0) stop 'Problem allocating t_dtmatJij%dtmat'
      ! else
      ! !allocate dtmat_xyz(lmmax,lmmax,3,irec) for irec_max=iemax_local
      ! allocate(t_dtmatJij_at%dtmat_xyz(t_inc%LMMAXD,t_inc%LMMAXD,3,t_mpi_c_grid%ntot2), STAT=ierr)
      ! if(ierr/=0) stop 'Problem allocating t_dtmatJij%dtmat for mpi'
      ! end if
      allocate(t_dtmatjij_at%dtmat_xyz(t_inc%lmmaxd,t_inc%lmmaxd,3,t_mpi_c_grid%ntot2),stat=i_stat)
      if (i_stat/=0) stop 'Problem allocating t_dtmatjij_at%dtmat_xyz'
      ! else
      ! allocate(t_tgmat%tmat(1,1,1), STAT=ierr)
      ! if(ierr/=0) stop 'Problem allocating dummy t_tgmat%tmat'
      ! end if

      t_dtmatjij_at%dtmat_xyz(:, :, :, :) = czero
    end if

  end subroutine init_t_dtmatjij_at

  !-------------------------------------------------------------------------------
  !> Summary: Store parameters needed in `t_imp`
  !> Author: 
  !> Category: initialization, KKRhost
  !> Deprecated: False 
  !> Store parameters needed in `t_imp`
  !-------------------------------------------------------------------------------
  subroutine init_params_t_imp(t_imp,ipand,natypd,irmd,irid,nfund,nspin,irmind,lmpotd)

    implicit none

    type (type_imp), intent (inout) :: t_imp

    integer, intent(in) :: irmd   !! Maximum number of radial points
    integer, intent(in) :: irid   !! Shape functions parameters in non-spherical part
    integer, intent(in) :: ipand  !! Number of panels in non-spherical part
    integer, intent(in) :: nfund  !! Shape functions parameters in non-spherical part
    integer, intent(in) :: nspin  !! Counter for spin directions
    integer, intent(in) :: natypd !! Number of kinds of atoms in unit cell
    integer, intent(in) :: irmind !! irmd - irnsd
    integer, intent(in) :: lmpotd !! (lpot+1)**2
    
    t_imp%ipand   = ipand
    t_imp%natypd  = natypd
    t_imp%irmd    = irmd
    t_imp%irid    = irid
    t_imp%nfund   = nfund
    t_imp%nspin   = nspin
    t_imp%irmind  = irmind
    t_imp%lmpotd  = lmpotd

  end subroutine init_params_t_imp

  !-------------------------------------------------------------------------------
  !> Summary: Subroutine to allocate and initialize arrays of `t_imp`
  !> Author: 
  !> Category: initialization, memory-management, KKRhost
  !> Deprecated: False 
  !> Subroutine to allocate and initialize arrays of `t_imp`
  !-------------------------------------------------------------------------------
  subroutine init_t_imp(t_inc, t_imp)

    implicit none

    type (type_inc), intent (in) :: t_inc
    type (type_imp), intent (inout) :: t_imp
    ! local
    integer :: irmd     !! Maximum number of radial points
    integer :: irid     !! Shape functions parameters in non-spherical part
    integer :: ipand    !! Number of panels in non-spherical part
    integer :: nfund    !! Shape functions parameters in non-spherical part
    integer :: nspin    !! Counter for spin directions
    integer :: natypd   !! Number of kinds of atoms in unit cell
    integer :: irmind   !! irmd - irnsd
    integer :: lmpotd   !! (lpot+1)**2
    integer :: natomimp !! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation

    integer :: i_stat
    ! so far only with SOC implemented
    if (.not. t_inc%newsosol) stop 'in init_t_imp: should only be called with NEWSOSOL'

    ! for convenience define this parameter locally
    natomimp  = t_imp%natomimp
    ipand     = t_imp%ipand
    natypd    = t_imp%natypd
    irmd      = t_imp%irmd
    irid      = t_imp%irid
    nfund     = t_imp%nfund
    nspin     = t_imp%nspin
    irmind    = t_imp%irmind
    lmpotd    = t_imp%lmpotd

    ! integer arrays
    if (.not. allocated(t_imp%ipanimp)) then
      allocate(t_imp%ipanimp(natomimp),stat=i_stat)
      if (i_stat/=0) stop 'Problem allocating t_imp%ipanimp'
      ! initialize with zeros
      t_imp%ipanimp = 0
    end if
    if (.not. allocated(t_imp%ircutimp)) then
      allocate(t_imp%ircutimp(0:ipand,natomimp),stat=i_stat)
      if (i_stat/=0) stop 'Problem allocating t_imp%ircutimp'
      t_imp%ircutimp = 0
    end if
    if (.not. allocated(t_imp%irminimp)) then
      allocate(t_imp%irminimp(natomimp),stat=i_stat)
      if (i_stat/=0) stop 'Problem allocating t_imp%irminimp'
      t_imp%irminimp = 0
    end if
    if (.not. allocated(t_imp%irwsimp)) then
      allocate(t_imp%irwsimp(natomimp),stat=i_stat)
      if (i_stat/=0) stop 'Problem allocating t_imp%irwsimp'
      t_imp%irwsimp = 0
    end if
    if (.not. allocated(t_imp%hostimp)) then
      allocate(t_imp%hostimp(natypd),stat=i_stat)
      if (i_stat/=0) stop 'Problem allocating t_imp%hostimp'
      t_imp%hostimp = 0
    end if
    if (.not. allocated(t_imp%atomimp)) then
      allocate(t_imp%atomimp(natomimp),stat=i_stat)
      if (i_stat/=0) stop 'Problem allocating t_imp%atomimp'
      t_imp%atomimp = 0
    end if

    ! real (kind=dp) arrays
    if (.not. allocated(t_imp%rimp)) then
      allocate(t_imp%rimp(irmd,natomimp),stat=i_stat)
      if (i_stat/=0) stop 'Problem allocating t_imp%rimp'
      t_imp%rimp = 0.0_dp
    end if
    if (.not. allocated(t_imp%zimp)) then
      allocate(t_imp%zimp(natomimp),stat=i_stat)
      if (i_stat/=0) stop 'Problem allocating t_imp%zimp'
      t_imp%zimp = 0.0_dp
    end if
    if (.not. allocated(t_imp%thetasimp)) then
      allocate(t_imp%thetasimp(irid,nfund,natomimp),stat=i_stat)
      if (i_stat/=0) stop 'Problem allocating t_imp%thetasimp'
      t_imp%thetasimp = 0.0_dp
    end if
    if (.not. allocated(t_imp%vispimp)) then
      allocate(t_imp%vispimp(irmd,natomimp*nspin),stat=i_stat)
      if (i_stat/=0) stop 'Problem allocating t_imp%vispimp'
      t_imp%vispimp = 0.0_dp
    end if
    if (.not. allocated(t_imp%vinsimp)) then
      allocate(t_imp%vinsimp(irmind:irmd,lmpotd,natomimp*nspin),stat=i_stat)
      if (i_stat/=0) stop 'Problem allocating t_imp%vispimp'
      t_imp%vinsimp = 0.0_dp
    end if
    if (.not. allocated(t_imp%rclsimp)) then
      allocate(t_imp%rclsimp(3,natomimp),stat=i_stat)
      if (i_stat/=0) stop 'Problem allocating t_imp%rclsimp'
      t_imp%rclsimp = 0.0_dp
    end if
    if (.not. allocated(t_imp%thetaimp)) then
      allocate(t_imp%thetaimp(natomimp),stat=i_stat)
      if (i_stat/=0) stop 'Problem allocating t_imp%thetaimp'
      t_imp%thetaimp = 0.0_dp
    end if
    if (.not. allocated(t_imp%phiimp)) then
      allocate(t_imp%phiimp(natomimp),stat=i_stat)
      if (i_stat/=0) stop 'Problem allocating t_imp%phiimp'
      t_imp%phiimp = 0.0_dp
    end if

    ! complex (kind=dp) arrays
    if (.not. allocated(t_imp%rllimp)) then
      allocate(t_imp%rllimp(t_inc%nsra*t_inc%lmmaxso,t_inc%lmmaxso,t_inc%irmdnew,natomimp),stat=i_stat)
      if (i_stat/=0) stop 'Problem allocating t_imp%rllimp'
      t_imp%rllimp = czero
    end if

  end subroutine init_t_imp


#ifdef CPP_MPI
  !-------------------------------------------------------------------------------
  !> Summary: Subroutine to broadcast `t_inc` and `t_tgmat` over mpi ranks
  !> Author: 
  !> Category: communication, memory-management, KKRhost
  !> Deprecated: False 
  !> Subroutine to broadcast `t_inc` and `t_tgmat` over mpi ranks
  !-------------------------------------------------------------------------------
  !> @note ruess: after myBcast_impcls from Pkkr_sidebranch2D_2014_12_16 by Bernd Zimmermann
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine bcast_t_inc_tgmat(t_inc, t_tgmat, t_cpa, master)

    use :: mpi
    implicit none

    type (type_inc), intent (inout) :: t_inc
    type (type_tgmatices), intent (inout) :: t_tgmat
    type (type_cpa), intent (inout) :: t_cpa

    integer, intent (in) :: master
    integer :: mympitype1 ! for parameter from t_inc
    integer, dimension(t_inc%nparams) :: blocklen1, etype1 ! for parameter from t_inc
    integer :: mympitype2 ! for logicals in t_tgmat
    integer, dimension(t_tgmat%nelements) :: blocklen2, etype2 ! for logicals in t_tgmat
    integer :: ierr
    integer (kind=mpi_address_kind) :: base
    integer (kind=mpi_address_kind), dimension(t_inc%nparams) :: disp1, disp2

    !--------------------------------------------------------------------------------
    ! broadcast parameters from t_inc
    !--------------------------------------------------------------------------------
    call mpi_get_address(t_inc%nparams, disp1(1), ierr)
    call mpi_get_address(t_inc%lmmaxd, disp1(2), ierr)
    call mpi_get_address(t_inc%nspin, disp1(3), ierr)
    call mpi_get_address(t_inc%ielast, disp1(4), ierr)
    call mpi_get_address(t_inc%nqdos, disp1(5), ierr)
    call mpi_get_address(t_inc%natyp, disp1(6), ierr)
    call mpi_get_address(t_inc%lmgf0d, disp1(7), ierr)
    call mpi_get_address(t_inc%nclsd, disp1(8), ierr)
    call mpi_get_address(t_inc%naclsmax, disp1(9), ierr)
    call mpi_get_address(t_inc%i_iteration, disp1(10), ierr)
    call mpi_get_address(t_inc%n_iteration, disp1(11), ierr)
    call mpi_get_address(t_inc%mit_bry, disp1(12), ierr)
    call mpi_get_address(t_inc%nshell0, disp1(13), ierr)
    call mpi_get_address(t_inc%nkmesh, disp1(14), ierr)
    call mpi_get_address(t_inc%newsosol, disp1(15), ierr)
    call mpi_get_address(t_inc%nosoc, disp1(16), ierr)
    call mpi_get_address(t_inc%deci_out, disp1(17), ierr)
    call mpi_get_address(t_inc%i_write, disp1(18), ierr)
    call mpi_get_address(t_inc%i_time, disp1(19), ierr)
    call mpi_get_address(t_inc%nsra, disp1(20), ierr)
    call mpi_get_address(t_inc%lmmaxso, disp1(21), ierr)
    call mpi_get_address(t_inc%irmdnew, disp1(22), ierr)
    call mpi_get_address(t_inc%kvrel, disp1(23), ierr)
    base = disp1(1)
    disp1 = disp1 - base

    blocklen1(1:23) = 1

    etype1(1:23) = mpi_integer
    etype1(15:17) = mpi_logical

    call mpi_type_create_struct(t_inc%nparams, blocklen1, disp1, etype1, mympitype1, ierr)
    if (ierr/=mpi_success) stop 'Problem in create_mpimask_t_inc'

    call mpi_type_commit(mympitype1, ierr)
    if (ierr/=mpi_success) stop 'error commiting create_mpimask_t_inc'

    call mpi_bcast(t_inc%nparams, 1, mympitype1, master, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop 'error brodcasting t_inc'

    call mpi_type_free(mympitype1, ierr)

    ! broadcast allocatable array kmesh(nkmesh)
    if (.not. allocated(t_inc%kmesh_ie)) allocate (t_inc%kmesh_ie(t_inc%ielast))
    if (.not. allocated(t_inc%kmesh)) allocate (t_inc%kmesh(t_inc%nkmesh))
    call mpi_bcast(t_inc%kmesh_ie, t_inc%ielast, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(t_inc%kmesh, t_inc%nkmesh, mpi_integer, master, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop 'error brodcasting t_inc%kmesh'
    !--------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------
    ! brodcast allocatable arrays from t_tgmat
    ! first broadcast logocal switches
    !--------------------------------------------------------------------------------
    call mpi_get_address(t_tgmat%nelements, disp2(1), ierr)
    call mpi_get_address(t_tgmat%tmat_to_file, disp2(2), ierr)
    call mpi_get_address(t_tgmat%gmat_to_file, disp2(3), ierr)
    call mpi_get_address(t_tgmat%gref_to_file, disp2(4), ierr)

    base = disp2(1)
    disp2 = disp2 - base

    blocklen2(1:4) = 1

    etype2(1) = mpi_integer
    etype2(2:4) = mpi_logical

    call mpi_type_create_struct(t_tgmat%nelements, blocklen2, disp2, etype2, mympitype2, ierr)
    if (ierr/=mpi_success) stop 'Problem in create_mpimask_tgmat_logicals'

    call mpi_type_commit(mympitype2, ierr)
    if (ierr/=mpi_success) stop 'error commiting create_mpimask_tgmat_logicals'

    call mpi_bcast(t_tgmat%nelements, 1, mympitype2, master, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop 'error brodcasting logicals from t_tgmat'

    call mpi_type_free(mympitype2, ierr)

    call mpi_bcast(t_cpa%dmatproj_to_file, 1, mpi_logical, master, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop 'error brodcasting logicals from t_cpa'

  end subroutine bcast_t_inc_tgmat
#endif


#ifdef CPP_MPI
  !-------------------------------------------------------------------------------
  !> Summary: Subroutine to broadcast Lloyd parameters over mpi ranks
  !> Author: 
  !> Category: communication, KKRhost
  !> Deprecated: False 
  !> Subroutine to broadcast Lloyd parameters over mpi ranks
  !-------------------------------------------------------------------------------
  subroutine bcast_t_lly_1(t_lloyd, master)

    use :: mpi
    implicit none

    type (type_lloyd), intent (inout) :: t_lloyd

    integer, intent (in) :: master
    integer :: mympitype1 ! for parameter from t_lloyd
    integer, dimension(t_lloyd%n1) :: blocklen1, etype1 ! for parameter from t_lloyd
    integer :: ierr
    integer (kind=mpi_address_kind) :: base
    integer (kind=mpi_address_kind), dimension(t_lloyd%n1) :: disp1

    call mpi_get_address(t_lloyd%n1, disp1(1), ierr)
    call mpi_get_address(t_lloyd%dtmat_to_file, disp1(2), ierr)
    call mpi_get_address(t_lloyd%tralpha_to_file, disp1(3), ierr)
    call mpi_get_address(t_lloyd%cdos_diff_lly_to_file, disp1(4), ierr)
    call mpi_get_address(t_lloyd%dgref_to_file, disp1(5), ierr)
    call mpi_get_address(t_lloyd%g0tr_to_file, disp1(6), ierr)

    base = disp1(1)
    disp1 = disp1 - base

    blocklen1(1:6) = 1

    etype1(1) = mpi_integer
    etype1(2:6) = mpi_logical

    call mpi_type_create_struct(t_lloyd%n1, blocklen1, disp1, etype1, mympitype1, ierr)
    if (ierr/=mpi_success) stop 'Problem in create_mpimask_tgmat_logicals'

    call mpi_type_commit(mympitype1, ierr)
    if (ierr/=mpi_success) stop 'error commiting create_mpimask_tgmat_logicals'

    call mpi_bcast(t_lloyd%n1, 1, mympitype1, master, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop 'error brodcasting logicals from t_tgmat'

    call mpi_type_free(mympitype1, ierr)

  end subroutine bcast_t_lly_1
#endif

  !-------------------------------------------------------------------------------
  !> Summary: Subroutine to allocate and initialize `t_lloyd`
  !> Author: 
  !> Category: initialization, memory-management, KKRhost
  !> Deprecated: False 
  !> Subroutine to allocate and initialize `t_lloyd`
  !-------------------------------------------------------------------------------
  subroutine init_tlloyd(t_inc, t_lloyd, t_mpi_c_grid)

    implicit none

    type (type_inc), intent (in) :: t_inc
    type (type_mpi_cartesian_grid_info), intent (in) :: t_mpi_c_grid
    type (type_lloyd), intent (inout) :: t_lloyd

    integer :: i_stat, nspin

    nspin = t_inc%nspin
    if (t_inc%newsosol .and. .not.t_inc%nosoc) nspin = 1  ! t_inc%NSPIN !1


    if (.not. allocated(t_lloyd%dtmat)) then
      if (.not. t_lloyd%dtmat_to_file) then
        ! if(nranks.eq.1) then
        ! !allocate dtmat(lmmax,lmmax,irec_max) for irec_max=ielast*nspin*natyp
        ! allocate(t_lloyd%dtmat(t_inc%LMMAXD,t_inc%LMMAXD,t_inc%IELAST*nspin*t_inc%NATYP), STAT=ierr)
        ! if(ierr/=0) stop 'Problem allocating t_tgmat%tmat'
        ! else
        ! !allocate dtmat(lmmax,lmmax,irec_max) for irec_max=iemax_local*nspin*natyp
        ! allocate(t_lloyd%dtmat(t_inc%LMMAXD,t_inc%LMMAXD,t_mpi_c_grid%ntot2*nspin*t_inc%NATYP), STAT=ierr)
        ! if(ierr/=0) stop 'Problem allocating t_tgmat%tmat for mpi'
        ! end if
        ! allocate dtmat(lmmax,lmmax,irec_max) for irec_max=iemax_local*nspin*natyp
        allocate(t_lloyd%dtmat(t_inc%lmmaxd,t_inc%lmmaxd,t_mpi_c_grid%ntot2*nspin*t_inc%natyp),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_lloyd%dtmat'
      else
        allocate(t_lloyd%dtmat(1,1,1),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_lloyd%dtmat'
      end if
      t_lloyd%dtmat(:, :, :) = czero
    end if

    if (.not. allocated(t_lloyd%tralpha)) then
      if (.not. t_lloyd%tralpha_to_file) then
        ! if(nranks.eq.1) then
        ! allocate(t_lloyd%tralpha(t_inc%IELAST*nspin*t_inc%NATYP), STAT=ierr)
        ! if(ierr/=0) stop 'Problem allocating t_lloyd%tralpha'
        ! else
        ! allocate(t_lloyd%tralpha(t_mpi_c_grid%ntot2*nspin*t_inc%NATYP), STAT=ierr)
        ! if(ierr/=0) stop 'Problem allocating t_lloyd%tralpha for mpi'
        ! end if
        allocate(t_lloyd%tralpha(t_mpi_c_grid%ntot2*nspin*t_inc%natyp),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_lloyd%tralpha'
      else
        allocate(t_lloyd%tralpha(1),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_lloyd%tralpha'
      end if
      t_lloyd%tralpha(:) = czero
    end if

    if (.not. allocated(t_lloyd%cdos)) then
      if (.not. t_lloyd%cdos_diff_lly_to_file) then
        ! if(nranks.eq.1) then
        ! allocate(t_lloyd%cdos(t_inc%IELAST,t_inc%NSPIN), STAT=ierr)
        ! if(ierr/=0) stop 'Problem allocating t_lloyd%cdos'
        ! else
        ! allocate(t_lloyd%cdos(t_inc%IELAST,t_inc%NSPIN), STAT=ierr)
        ! if(ierr/=0) stop 'Problem allocating t_lloyd%cdos for mpi'
        ! end if
        allocate(t_lloyd%cdos(t_inc%ielast,t_inc%nspin),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_lloyd%cdos'
      else
        allocate(t_lloyd%cdos(1,1),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_lloyd%cdos'
      end if
      t_lloyd%cdos(:, :) = czero
    end if

    if (.not. allocated(t_lloyd%dgref)) then
      if (.not. t_lloyd%dgref_to_file) then
        ! if(nranks.eq.1) then
        ! allocate(t_lloyd%dgref(t_inc%NACLSMAX*t_inc%LMGF0D,t_inc%LMGF0D,t_inc%NCLSD,t_inc%IELAST), STAT=ierr)
        ! if(ierr/=0) stop 'Problem allocating t_lloyd%dgref'
        ! else
        ! allocate(t_lloyd%dgref(t_inc%NACLSMAX*t_inc%LMGF0D,t_inc%LMGF0D,t_inc%NCLSD,t_mpi_c_grid%ntot2), STAT=ierr)
        ! if(ierr/=0) stop 'Problem allocating t_lloyd%dgref for mpi'
        ! end if
        allocate(t_lloyd%dgref(t_inc%naclsmax*t_inc%lmgf0d,t_inc%lmgf0d,t_inc%nclsd,t_mpi_c_grid%ntot2),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_lloyd%dgref'
      else
        allocate(t_lloyd%dgref(1,1,1,1),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_lloyd%dgref'
      end if
      t_lloyd%dgref(:, :, :, :) = czero
    end if

    if (.not. allocated(t_lloyd%g0tr)) then
      if (.not. t_lloyd%g0tr_to_file) then
        ! if(nranks.eq.1) then
        ! allocate(t_lloyd%g0tr(t_inc%IELAST), STAT=ierr)
        ! if(ierr/=0) stop 'Problem allocating t_lloyd%g0tr'
        ! else
        ! allocate(t_lloyd%g0tr(t_mpi_c_grid%ntot2), STAT=ierr)
        ! if(ierr/=0) stop 'Problem allocating t_lloyd%g0tr for mpi'
        ! end if
        allocate(t_lloyd%g0tr(t_mpi_c_grid%ntot2),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_lloyd%g0tr'
      else
        allocate(t_lloyd%g0tr(1),stat=i_stat)
        if (i_stat/=0) stop 'Problem allocating t_lloyd%g0tr'
      end if
      t_lloyd%g0tr(:) = czero
    end if

  end subroutine init_tlloyd

#ifdef CPP_MPI
  !-------------------------------------------------------------------------------
  !> Summary: Subroutine to store MPI rank for of 2 level parallelization
  !> Author: 
  !> Category: initialization, communication, KKRhost
  !> Deprecated: False 
  !> Subroutine to store MPI rank for of 2 level parallelization
  !-------------------------------------------------------------------------------
  subroutine save_t_mpi_c_grid(t_mpi_c_grid,subarr_dim,mympi_comm_ie,mympi_comm_at, &
    myrank_ie,myrank_at,myrank_atcomm,nranks_ie,nranks_at,nranks_atcomm)

    use :: mpi
    implicit none
    type (type_mpi_cartesian_grid_info), intent (inout) :: t_mpi_c_grid
    integer, dimension(2), intent (in) :: subarr_dim
    integer, intent (in) :: mympi_comm_ie, mympi_comm_at, myrank_ie, myrank_at, nranks_ie, nranks_at, nranks_atcomm, myrank_atcomm

    t_mpi_c_grid%dims           = subarr_dim
    t_mpi_c_grid%mympi_comm_ie  = mympi_comm_ie
    t_mpi_c_grid%mympi_comm_at  = mympi_comm_at
    t_mpi_c_grid%myrank_ie      = myrank_ie
    t_mpi_c_grid%myrank_at      = myrank_at
    t_mpi_c_grid%myrank_atcomm  = myrank_atcomm
    t_mpi_c_grid%nranks_ie      = nranks_ie
    t_mpi_c_grid%nranks_at      = nranks_at
    t_mpi_c_grid%nranks_atcomm  = nranks_atcomm

  end subroutine save_t_mpi_c_grid
#endif


#ifdef CPP_MPI
  !-------------------------------------------------------------------------------
  !> Summary: Subroutine to extract number of elements and offsets from `t_mpi_c_grid`
  !> Author: 
  !> Category: initialization, communication, KKRhost
  !> Deprecated: False 
  !> Subroutine to extract number of elements and offsets from `t_mpi_c_grid`
  !-------------------------------------------------------------------------------
  subroutine get_ntot_pt_ioff_pt_2d(t_mpi_c_grid, ntot_all, ioff_all)

    use :: mpi
    implicit none
    type (type_mpi_cartesian_grid_info), intent (in) :: t_mpi_c_grid
    integer, dimension(t_mpi_c_grid%nranks_ie*t_mpi_c_grid%nranks_at), intent (out) :: ntot_all
    integer, dimension(t_mpi_c_grid%nranks_ie*t_mpi_c_grid%nranks_at), intent (out) :: ioff_all

    integer, dimension(t_mpi_c_grid%nranks_ie) :: ntot_pt1, ioff_pt1
    integer, dimension(t_mpi_c_grid%nranks_at) :: ntot_pt2, ioff_pt2
    integer :: n1, n2, i1, i2, i3

    ntot_pt1 = t_mpi_c_grid%ntot_pt1
    ioff_pt1 = t_mpi_c_grid%ioff_pt1
    ntot_pt2 = t_mpi_c_grid%ntot_pt2
    ioff_pt2 = t_mpi_c_grid%ioff_pt2
    n1 = t_mpi_c_grid%nranks_ie
    n2 = t_mpi_c_grid%nranks_at


    do i1 = 1, n1
      do i2 = 1, n2
        i3 = i2 + n2*(i1-1)
        ntot_all(i3) = ntot_pt1(i1)*ntot_pt2(i2)
        if (i3==1) then
          ioff_all(i3) = 0
        else
          ioff_all(i3) = ioff_all(i3-1) + ntot_all(i3-1)
        end if
      end do
    end do

  end subroutine get_ntot_pt_ioff_pt_2d
#endif


#ifdef CPP_MPI
  !-------------------------------------------------------------------------------
  !> Summary: Subroutine to communicate arrays in `t_lloyd` over ranks
  !> Author: 
  !> Category: memory-management, communication, KKRhost
  !> Deprecated: False 
  !> Subroutine to communicate arrays in `t_lloyd` over ranks
  !-------------------------------------------------------------------------------
  subroutine gather_lly_dtmat(t_mpi_c_grid, t_lloyd, lmmaxd, mympi_comm)

    use :: mpi
    implicit none

    type (type_mpi_cartesian_grid_info), intent (in) :: t_mpi_c_grid
    type (type_lloyd), intent (inout) :: t_lloyd
    integer, intent (in) :: lmmaxd, mympi_comm

    complex (kind=dp), dimension(:,:,:), allocatable :: work_lly
    integer :: ierr,i_stat,i_all, iwork, nspin

    nspin = t_inc%nspin
    if (t_inc%newsosol .and. .not.t_inc%nosoc) nspin = 1


    ! communicate dtmatll
    iwork = lmmaxd*lmmaxd*t_mpi_c_grid%ntot2*nspin*t_inc%natyp
    if (iwork/=product(shape(t_lloyd%dtmat))) stop '[gather_lly_dtmat]: Error shape mismatch in gather_dtmat_lly'
 
    allocate(work_lly(lmmaxd,lmmaxd,t_mpi_c_grid%ntot2*nspin*t_inc%natyp),stat=i_stat)
    if (i_stat/=0) stop 'Problem allocating work_lly'
 
    call mpi_allreduce(t_lloyd%dtmat, work_lly, iwork, mpi_double_complex, mpi_sum, mympi_comm, ierr)
    call zcopy(iwork, work_lly, 1, t_lloyd%dtmat, 1)

    deallocate (work_lly, stat=i_stat)
    if (i_stat/=0) stop 'Problem deallocating work_lly'

    ! communicate tralpha
    iwork = t_mpi_c_grid%ntot2*nspin*t_inc%natyp
    if (iwork/=product(shape(t_lloyd%tralpha))) stop '[gather_lly_dtmat]: Error shape mismatch in gather_dtmat_lly'

    allocate(work_lly(iwork,1,1),stat=i_stat)
    if (i_stat/=0) stop 'Problem allocating work_lly'

    call mpi_allreduce(t_lloyd%tralpha, work_lly(:,1,1), iwork, mpi_double_complex, mpi_sum, mympi_comm, ierr)
    call zcopy(iwork, work_lly(:,1,1), 1, t_lloyd%tralpha, 1)

    deallocate (work_lly, stat=i_stat)
    if (i_stat/=0) stop 'Problem deallocating work_lly'
  end subroutine gather_lly_dtmat
#endif


#ifdef CPP_MPI
  !-------------------------------------------------------------------------------
  !> Summary: Subroutine to communicate single site t-matrix over ranks
  !> Author: 
  !> Category: memory-management, communication, single-site, KKRhost
  !> Deprecated: False 
  !> Subroutine to communicate single site t-matrix over ranks
  !-------------------------------------------------------------------------------
  subroutine gather_tmat(t_inc, t_tgmat, t_mpi_c_grid, ntot_pt, ioff_pt, mytot, mympi_comm, nranks)

    use :: mpi
    implicit none

    type (type_inc), intent (in) :: t_inc
    type (type_mpi_cartesian_grid_info), intent (in) :: t_mpi_c_grid
    type (type_tgmatices), intent (inout) :: t_tgmat
    integer, intent (in) :: nranks
    integer, intent (in) :: mytot, mympi_comm
    integer, dimension(0:nranks-1), intent (in) :: ntot_pt, ioff_pt

    integer :: ihelp
    integer :: nspin
    integer, dimension(nranks) :: recvcounts, displs
    complex (kind=dp), dimension(:,:,:), allocatable :: work
    integer :: ierr, idim, i_stat,i_all


    ! Gather tmat so that all processors the full matrix for their part of the energy contour
    if (t_mpi_c_grid%nranks_ie>1) then
      nspin = t_inc%nspin
      if (t_inc%newsosol .and. .not.t_inc%nosoc) nspin = 1

      ihelp = t_inc%lmmaxd**2*t_mpi_c_grid%ntot2*nspin ! *t_inc%NATYP/mytot
      recvcounts = ntot_pt*ihelp
      displs = ioff_pt*ihelp

      allocate(work(t_inc%lmmaxd,t_inc%lmmaxd,t_mpi_c_grid%ntot2*nspin*t_inc%natyp),stat=i_stat)
      if (i_stat/=0) stop 'Problem allocating work'
    
      call mpi_allgatherv(t_tgmat%tmat, mytot*ihelp, mpi_double_complex, work, recvcounts, displs, mpi_double_complex, mympi_comm, ierr)
      idim = t_inc%lmmaxd**2*t_mpi_c_grid%ntot2*nspin*t_inc%natyp
      call zcopy(idim, work, 1, t_tgmat%tmat, 1)

      deallocate (work, stat=i_stat)
      if (i_stat/=0) stop 'Problem deallocating work'
    end if

  end subroutine gather_tmat
#endif


#ifdef CPP_MPI
  ! subroutine gather_gref(t_inc, t_tgmat, t_mpi_c_grid, ntot_pT, ioff_pT, mytot, mympi_comm)

  ! use mpi
  ! use mod_mympi,   only: myrank, master
  ! implicit none

  ! type(type_inc), intent(in) :: t_inc
  ! type(type_mpi_cartesian_grid_info), intent(in) :: t_mpi_c_grid
  ! type(type_tgmatices), intent(inout) :: t_tgmat
  ! integer, intent(in) :: ntot_pT(0:t_mpi_c_grid%nranks_ie-1), ioff_pT(0:t_mpi_c_grid%nranks_ie-1), mytot, mympi_comm

  ! integer :: ihelp
  ! integer :: ierr

  ! !Gather gref so that all processors have the full matrix for their part of the energy countour
  ! if(t_mpi_c_grid%dims(1)>1) then
  ! ihelp      = t_inc%NACLSD*t_inc%LMGF0D*t_inc%LMGF0D*t_inc%NCLSD
  ! call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  ! call MPI_Bcast( t_tgmat%gref, mytot*ihelp, MPI_DOUBLE_COMPLEX, 0 , &
  ! & mympi_comm, ierr )
  ! end if

  ! end subroutine gather_gref
#endif


#ifdef CPP_MPI
  !-------------------------------------------------------------------------------
  !> Summary: Subroutine to communicate structural green function over ranks
  !> Author: 
  !> Category: communication, structural-greensfunction, KKRhost
  !> Deprecated: False 
  !> Subroutine to communicate structural green function over ranks
  !-------------------------------------------------------------------------------
  subroutine gather_gmat(t_inc, t_tgmat, ntot_pt, ioff_pt, mytot, nranks)

    use :: mpi
    implicit none

    type (type_inc), intent (in) :: t_inc
    type (type_tgmatices), intent (inout) :: t_tgmat
    integer, intent (in) :: nranks
    integer, intent (in) :: mytot
    integer, dimension(0:nranks-1), intent (in) :: ntot_pt, ioff_pt

    integer :: ihelp
    integer :: recvcounts(0:nranks-1), displs(0:nranks-1)
    integer :: ierr, nspin

    nspin = t_inc%nspin
    if (t_inc%newsosol .and. .not.t_inc%nosoc) nspin = 1

    ! Gather gmat so that all processors have the full matrix
    ihelp = t_inc%lmmaxd*t_inc%lmmaxd*t_inc%nqdos*nspin ! *t_inc%IELAST*t_inc%NSPIN*t_inc%NATYP
    if (t_mpi_c_grid%dims(1)>1) then
      recvcounts = ntot_pt*ihelp
      displs = ioff_pt*ihelp
      call mpi_allgatherv(t_tgmat%gmat, mytot*ihelp, mpi_double_complex, t_tgmat%gmat, recvcounts, displs, mpi_double_complex, mpi_comm_world, ierr)
    end if

  end subroutine gather_gmat
#endif


#ifdef CPP_MPI
  !-------------------------------------------------------------------------------
  !> Summary: Subroutine to communicate scalars of `t_imp` over ranks
  !> Author: 
  !> Category: communication, KKRhost
  !> Deprecated: False 
  !> Subroutine to communicate scalars of `t_imp` over ranks
  !-------------------------------------------------------------------------------
  subroutine bcast_t_imp_scalars(t_imp, master)

    use :: mpi
    implicit none

    type (type_imp), intent (inout) :: t_imp

    integer, intent (in) :: master
    integer :: mympitype1 ! for scalars from t_imp
    integer, dimension(t_imp%n1) :: blocklen1, etype1 ! for scalars from t_imp
    integer :: ierr
    integer (kind=mpi_address_kind) :: base
    integer (kind=mpi_address_kind), dimension(t_imp%n1) :: disp1

    ! scalars
    call mpi_get_address(t_imp%n1, disp1(1), ierr)
    call mpi_get_address(t_imp%n2, disp1(2), ierr)
    call mpi_get_address(t_imp%natomimp, disp1(3), ierr)
    call mpi_get_address(t_imp%ihost, disp1(4), ierr)
    call mpi_get_address(t_imp%ipand, disp1(5), ierr)
    call mpi_get_address(t_imp%natypd, disp1(6), ierr)
    call mpi_get_address(t_imp%irmd, disp1(7), ierr)
    call mpi_get_address(t_imp%irid, disp1(8), ierr)
    call mpi_get_address(t_imp%nfund, disp1(9), ierr)
    call mpi_get_address(t_imp%nspin, disp1(10), ierr)
    call mpi_get_address(t_imp%irmind, disp1(11), ierr)
    call mpi_get_address(t_imp%lmpotd, disp1(12), ierr)

    base = disp1(1)
    disp1 = disp1 - base

    blocklen1(:) = 1

    etype1(:) = mpi_integer

    call mpi_type_create_struct(t_imp%n1, blocklen1, disp1, etype1, mympitype1, ierr)
    if (ierr/=mpi_success) stop 'Problem in create_mpimask_tgmat_logicals'
    call mpi_type_commit(mympitype1, ierr)
    if (ierr/=mpi_success) stop 'error commiting create_mpimask_tgmat_logicals'
    call mpi_bcast(t_imp%n1, 1, mympitype1, master, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop 'error brodcasting logicals from t_tgmat'
    call mpi_type_free(mympitype1, ierr)

  end subroutine bcast_t_imp_scalars

  !-------------------------------------------------------------------------------
  !> Summary: Subroutine to communicate arrays of t_imp over ranks
  !> Author: 
  !> Category: communication, KKRhost
  !> Deprecated: False 
  !> Subroutine to communicate arrays of t_imp over ranks
  !-------------------------------------------------------------------------------
  subroutine bcast_t_imp_arrays(t_imp, t_inc, master)

    use :: mpi
    implicit none

    type (type_imp), intent (inout) :: t_imp
    type (type_inc), intent (in) :: t_inc

    integer, intent (in) :: master
    integer :: mympitype2 ! for arrays from t_imp
    integer, dimension(t_imp%n2) :: blocklen2, etype2 ! for arrays from t_imp
    integer :: ierr
    integer (kind=mpi_address_kind) :: base
    integer (kind=mpi_address_kind), dimension(t_imp%n2) :: disp2

    ! arrays
    call mpi_get_address(t_imp%n2, disp2(1), ierr)
    call mpi_get_address(t_imp%ipanimp, disp2(2), ierr)
    call mpi_get_address(t_imp%ircutimp, disp2(3), ierr)
    call mpi_get_address(t_imp%irminimp, disp2(4), ierr)
    call mpi_get_address(t_imp%irwsimp, disp2(5), ierr)
    call mpi_get_address(t_imp%hostimp, disp2(6), ierr)
    call mpi_get_address(t_imp%rimp, disp2(7), ierr)
    call mpi_get_address(t_imp%zimp, disp2(8), ierr)
    call mpi_get_address(t_imp%thetasimp, disp2(9), ierr)
    call mpi_get_address(t_imp%vispimp, disp2(10), ierr)
    call mpi_get_address(t_imp%vinsimp, disp2(11), ierr)
    call mpi_get_address(t_imp%rclsimp, disp2(12), ierr)
    call mpi_get_address(t_imp%rllimp, disp2(13), ierr)
    call mpi_get_address(t_imp%atomimp, disp2(14), ierr)
    call mpi_get_address(t_imp%thetaimp, disp2(15), ierr)
    call mpi_get_address(t_imp%phiimp, disp2(16), ierr)


    base = disp2(1)
    disp2 = disp2 - base

    blocklen2(1) = 1
    blocklen2(2) = t_imp%natomimp
    blocklen2(3) = (1+t_imp%ipand)*t_imp%natomimp
    blocklen2(4) = t_imp%natomimp
    blocklen2(5) = t_imp%natomimp
    blocklen2(6) = t_imp%natypd
    blocklen2(7) = t_imp%irmd*t_imp%natomimp
    blocklen2(8) = t_imp%natomimp
    blocklen2(9) = t_imp%irid*t_imp%nfund*t_imp%natomimp
    blocklen2(10) = t_imp%irmd*t_imp%natomimp*t_imp%nspin
    blocklen2(11) = (t_imp%irmd-t_imp%irmind+1)*t_imp%lmpotd*t_imp%natomimp*t_imp%nspin
    blocklen2(12) = 3*t_imp%natomimp
    blocklen2(13) = t_inc%nsra*t_inc%lmmaxso*t_inc%lmmaxso*t_inc%irmdnew*t_imp%natomimp
    blocklen2(14) = t_imp%natomimp
    blocklen2(15) = t_imp%natomimp
    blocklen2(16) = t_imp%natomimp

    etype2(1:6) = mpi_integer
    etype2(7:12) = mpi_double_precision
    etype2(13) = mpi_double_complex
    etype2(14) = mpi_integer
    etype2(15:16) = mpi_double_precision

    call mpi_type_create_struct(t_imp%n2, blocklen2, disp2, etype2, mympitype2, ierr)
    if (ierr/=mpi_success) stop 'Problem in create_mpimask_tgmat_logicals'
    call mpi_type_commit(mympitype2, ierr)
    if (ierr/=mpi_success) stop 'error commiting create_mpimask_tgmat_logicals'
    call mpi_bcast(t_imp%n2, 1, mympitype2, master, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop 'error brodcasting logicals from t_tgmat'
    call mpi_type_free(mympitype2, ierr)

  end subroutine bcast_t_imp_arrays
#endif


end module mod_types
