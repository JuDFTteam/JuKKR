!> multiple scattering k-loop and symmetrisation
module kloopz1_mod
  use SolverOptions_mod, only:
  implicit none
  private
  public :: kloopz1_new

  double complex, parameter :: cone=(1.d0,0.d0), zero=(0.d0,0.d0)
  
  CONTAINS

    subroutine KLOOPZ1_new(GMATN, solv, kkr_op, precond, ALAT, NOFKS, VOLBZ, BZKP, VOLCUB, RR, GINP_LOCAL, &
                           NSYMAT,DSYMLL, TMATLL, lmmaxd, nrd, trunc2atom_index, communicator, iguess_data)

! **********************************************************************

! Only part of arrays for corresponding spin direction is passed
! (GMATN, TSST_LOCAL, DTDE_LOCAL, LLY_GRDT, TR_ALPH, GMATXIJ)
!
! NOFKS .. number of k-points, integer
! VOLBZ .. Brillouin zone volume, double
! BZKP ... k-points of used k-mesh ... dimension (3, KPOIBZ)
! VOLCUB . array of Brillouin zone integration weights for each k-point ... dimension (KPOIBZ)

! GINP_LOCAL ... reference Green's function
! DGINP ...      derivative of reference Green's function
! TSST_LOCAL ..  t-matrix

    use kkrmat_new_mod, only: KKRMAT01_new
    use TEST_lcutoff_mod, only: cutoffmode
    use InitialGuess_mod, only: InitialGuess
    use ClusterInfo_mod, only: ClusterInfo
    use TFQMRSolver_mod, only: TFQMRSolver
    use BCPOperator_mod, only: BCPOperator
    use KKROperator_mod, only: KKROperator
    use MultScatData_mod, only: MultScatData
    use jij_calc_mod, only: global_jij_data, symjij
    integer, parameter :: NSYMAXD = 48

    class (TFQMRSolver), intent(inout) :: solv
    class (KKROperator), intent(inout) :: kkr_op
    class (BCPOperator), intent(inout) :: precond

    integer, intent(in) :: lmmaxd
    integer, intent(in) :: nrd
    !> mapping trunc. index -> atom index
    integer, intent(in) :: trunc2atom_index(:)
    integer, intent(in) :: communicator
    type(InitialGuess), intent(inout) :: iguess_data
    double precision, intent(in) :: ALAT
    double precision, intent(in) :: VOLBZ
    integer, intent(in) :: NOFKS
    integer, intent(in) :: NSYMAT
    double complex, intent(in) :: DSYMLL(LMMAXD,LMMAXD,NSYMAXD)
    double complex, intent(inout) :: GMATN(:,:,:)
    double complex, intent(inout) :: GINP_LOCAL(:,:,:,:)
    double complex, intent(inout) :: TMATLL(:,:,:)
    double precision, intent(in) :: RR(3,0:NRD)
    double precision, intent(in) :: BZKP(:,:) ! dim (3,kpoibz)
    double precision, intent(in) :: VOLCUB(:) ! dim kpoibz
    
    double complex :: TAUVBZ
    double precision :: RFCTOR
    integer :: INFO    ! for LAPACK calls
    integer :: symmetry_index
    integer :: ispin
    double complex :: GLL(LMMAXD,LMMAXD)
    double complex, allocatable :: GS(:,:,:,:)
    !     effective (site-dependent) Delta_t^(-1) matrix
    double complex, allocatable :: MSSQ (:,:,:)
    integer :: IPVT(LMMAXD) ! work array for LAPACK
    double complex :: work_array(LMMAXD,LMMAXD) ! work array for LAPACK ZGETRI
    double complex ::        TPG(LMMAXD,LMMAXD)
    double complex ::         XC(LMMAXD,LMMAXD)      ! to store temporary matrix-matrix mult. result
    integer :: num_local_atoms
    integer :: ilocal
    integer :: memory_stat

    type(MultScatData), pointer :: ms
    type(ClusterInfo), pointer :: cluster_info

    ms => kkr_op%get_ms_workspace()
    cluster_info => ms%cluster_info

    num_local_atoms = size(ms%atom_indices)

! -------------------------------------------------------------------
! Allocate Arrays
! -------------------------------------------------------------------
    allocate(GS(LMMAXD,LMMAXD,NSYMAXD,num_local_atoms), MSSQ(LMMAXD,LMMAXD,num_local_atoms), stat=memory_stat)

    if (memory_stat /= 0) then
      write(*,*) "KLOOPZ1: FATAL Error, failure to allocate memory, probably out of memory."
      stop
    endif

!     RFCTOR=A/(2*PI) conversion factor to p.u.
    RFCTOR = ALAT/(8.D0*ATAN(1.0D0))           ! = ALAT/(2*PI)

    do ilocal = 1, num_local_atoms

      MSSQ(1:LMMAXD,1:LMMAXD,ilocal) = TMATLL(1:LMMAXD,1:LMMAXD,ms%atom_indices(ilocal))

  ! ---> inversion

  !     The (local) Delta_t matrix is inverted and stored in MSSQ

      call ZGETRF(LMMAXD,LMMAXD,MSSQ(:,:,ilocal),LMMAXD,IPVT,INFO)
      call ZGETRI(LMMAXD,MSSQ(:,:,ilocal),LMMAXD,IPVT,work_array, LMMAXD*LMMAXD,INFO)

    enddo ! ilocal

!=======================================================================
!     Note: the actual k-loop is in kkrmat01 (it is not parallelized)
!     The integration over k is also performed in kkrmat01

    TAUVBZ = 1.D0/VOLBZ
    ! cutoffmode:
    ! 0 no cutoff
    ! 1 T-matrix cutoff, 2 full matrix cutoff (not supported anymore)
    ! 3 T-matrix cutoff with new solver
    ! 4 T-matrix cutoff with direct solver
    if (cutoffmode > 2 .or. cutoffmode == 0) then
      call KKRMAT01_new(solv, kkr_op, precond, BZKP,NOFKS,GS,VOLCUB,TMATLL, ALAT, NSYMAT, RR, GINP_LOCAL, lmmaxd, trunc2atom_index, communicator, iguess_data)
    else
      write(*,*) "0 < cutoffmode < 3 not supported."
      STOP
    endif
!-------------------------------------------------------- SYMMETRISE GLL


!      kkrmat01 returns GS (local) which contains NSYMAT -copies- (!)
!      (see 3rd index) of the scattering path operator
!      (already integrated over the irreducible wedge in k-space)
!      scattering path operator: ((Delta_T)^-1 - G_ref)^-1

!      All the symmetry operations are applied on GS and summed over
!     - the result is stored in GLL
!      Note: the symmetry operations apply on the (LL')-space

!------------------------------------------------------------------------------
    do ilocal = 1, num_local_atoms
!------------------------------------------------------------------------------

      do symmetry_index = 1,NSYMAT

    ! --->    GLL = sum(i=1,iumax)(tauvbz * DLL(i) * GS * DLL(i)^H)

        if (symmetry_index == 1) then

        ! --->    ull(1) is equal to unity matrix

          call ZCOPY(LMMAXD*LMMAXD,GS(1,1,1,ilocal),1,GLL,1)
          call ZSCAL(LMMAXD*LMMAXD,TAUVBZ,GLL,1)

        else

        ! --->      tpg = tauvbz * DLL * GS

          call ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,TAUVBZ, DSYMLL(1,1,symmetry_index),LMMAXD,GS(1,1,symmetry_index,ilocal),LMMAXD, zero,TPG,LMMAXD)

        ! --->    GLL = GLL + TPG * DLL(i)^H
        !                           C  ! dsymll might be complex in REL case

          call ZGEMM('N','C',LMMAXD,LMMAXD,LMMAXD,CONE,TPG,LMMAXD, DSYMLL(1,1,symmetry_index),LMMAXD,CONE,GLL,LMMAXD)
        endif ! symmetry_index == 1

      enddo ! symmetry_index
  !-------------------------------------------------------- IU = 1,NSYMAT


  ! --->  XC = Delta_t^(-1) * GLL

      call ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,MSSQ(:,:,ilocal), LMMAXD,GLL,LMMAXD,zero,XC,LMMAXD)

  !       GLL is overwritten with the following expression:
  !       (for the in configuration space diagonal structural Green's
  !       function of the REAL system - already integrated over k)

  ! --->  GLL = - Delta_t^(-1) - Delta_t^(-1) * GLL * Delta_t^(-1)

  !       copy overwrite GLL with content of MSSQ
      call ZCOPY(LMMAXD**2,MSSQ(:,:,ilocal),1,GLL,1)


  !       GLL = (-1) * XC                     *    MSSQ     + (-1) * GLL
  !                    |                            |                 |
  !            Delta_t^-1 * scat. path op.      Delta_t^-1    Delta_t^-1


!       call ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,-CONE,XC,LMMAXD, MSSQ(:,:,ilocal),LMMAXD,-CONE,GLL,LMMAXD)

      GLL = - XC - GLL

  ! --->  GMATN = GMATLL = GLL/RFCTOR...............rescaled and copied into output array

      GMATN(1:LMMAXD,1:LMMAXD,ilocal) = GLL(1:LMMAXD,1:LMMAXD)/RFCTOR

!------------------------------------------------------------------------------
    enddo ! ilocal
!------------------------------------------------------------------------------

    if (global_jij_data%do_jij_calculation) then

      ispin = global_jij_data%active_spin

      call SYMJIJ(ALAT,TAUVBZ, NSYMAT,DSYMLL, global_jij_data%NXIJ,global_jij_data%IXCP, &
                  TMATLL,MSSQ, global_jij_data%GSXIJ, global_jij_data%GMATXIJ(:,:,:,ispin), &  ! Result
                  cluster_info%naez_trc, lmmaxd, global_jij_data%nxijd)
    endif ! jij

! -------------------------------------------------------------------
! Deallocate Arrays
! -------------------------------------------------------------------

    deallocate(GS, MSSQ)

  endsubroutine KLOOPZ1_new

endmodule kloopz1_mod
