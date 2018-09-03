! -------------------------------------------------------------------------------
! MODULE: global_variables
! > @brief Module containing the variables which were previously in the inc.p
! > file
! > @author Jonathan Chico
! -------------------------------------------------------------------------------
module global_variables

  implicit none

  integer :: irid                  ! < Shape functions parameters in
                                   ! non-spherical part
  integer :: krel = 0              ! < Switch for
                                   ! non-relativistic/relativistic (0/1)
                                   ! program. Attention: several other
                                   ! parameters depend explicitly on KREL,
                                   ! they are set automatically Used for Dirac
                                   ! solver in ASA
  integer :: nfund                 ! < Shape functions parameters in
                                   ! non-spherical part
  integer :: ipand                 ! < Number of panels in non-spherical part
  integer :: ngshd                 ! < Shape functions parameters in
                                   ! non-spherical part
  integer :: ncleb                 ! < Number of Clebsch-Gordon coefficients
  integer :: knoco = 0             ! < (0/1) Collinear/Non-collinear magnetism
                                   ! (even in non-relativistic non-spin-orbit
                                   ! case)
  integer :: iemxd = 101           ! < Dimension for energy-dependent arrays
  integer :: irnsd = 890
  integer :: nmaxd = 2000000       ! < Paremeters for the Ewald summations
  integer :: ishld = 200000        ! < Paremeters for the Ewald summations
  integer :: naclsd                ! < Maximum number of atoms in a TB-cluster
  integer :: nspotd                ! < Number of potentials for storing
                                   ! non-sph. potentials
  integer :: ntperd                ! < Parameter in broyden subroutines
  integer :: ntrefd = 0            ! < parameter in broyden subroutine MUST BE
                                   ! 0 for the host program
  integer :: nsheld = 1000         ! < Number of blocks of the GF matrix that
                                   ! need to be calculated (NATYPD +
                                   ! off-diagonals in case of impurity)
  integer :: ncelld                ! < Number of cells (shapes) in
                                   ! non-spherical part
  integer :: nspind                ! < KREL+(1-KREL)*(NSPIN+1)
  integer :: knosph = 1            ! < switch for spherical/non-spherical
                                   ! (0/1) program.
  integer :: korbit = 0            ! < Spin-orbit/non-spin-orbit (1/0) added
                                   ! to the Schroedinger or SRA equations.
                                   ! Works with FP. KREL and KORBIT cannot be
                                   ! both non-zero.
  integer :: kpoibz = 250000       ! < Number of reciprocal space vectors
  integer :: wlength               ! < Word length for direct access files,
                                   ! compiler dependent ifort/others (1/4)
  integer :: nprincd = 1           ! < Number of principle layers, set to a
                                   ! number >= NRPINC in output of main0
  integer :: nlayerd               ! < Number of principal layers
                                   ! (NAEZD/NPRINCD) used in the inversion
                                   ! routines (independent on NATYPD)
  integer :: natomimpd = 150       ! < Size of the cluster for
                                   ! impurity-calculation output of GF should
                                   ! be 1, if you don't do such a calculation
  logical :: lnc                   ! < Coupled equations in two spins
                                   ! (switches true if KREL=1 or KORBIT=1 or
                                   ! KNOCO=1)

  integer :: lmaxd = 3             ! < lmax cutoff
  integer :: lmmaxd                ! < (KREL+KORBIT+1)*(LMAXD+1)^2
  integer :: lmgf0d                ! < (lmaxd+1)^2
  ! *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     * 
  ! *          function, set up in the spin-independent non-relativstic * 
  ! *          (l,m_l)-representation                                   * 
  ! *                                                                   * 
  ! ********************************************************************* 
  integer :: alm                   ! < naez*lmmaxd
  integer :: almgf0                ! < naezd*lmgf0d (<=alm)
  integer :: ndim_slabinv          ! < nprincd*lmmaxd
  integer :: nembd = 20            ! <
  integer :: nembd1 = 21           ! < NEBMD+1
  integer :: nembd2 = 21           ! < NEBMD+NAEZ
  integer :: nrd = 20000
  integer :: lm2d
  integer :: nclsd
  integer :: mmaxd
  integer :: npotd
  integer :: lmxspd
  integer :: lassld
  integer :: irmind
  integer :: nofgij
  integer :: nspindd
  integer :: nsatypd
  integer :: nrefd = 1
  integer :: irmd = 900
  integer :: naezd = 1
  integer :: natypd = 1
  integer :: lmpotd
  integer :: ntotd
  integer :: nrmaxd
  integer :: lmmaxso
  integer :: lpotd
  integer :: nchebd

  !> maximal number of different k-meshes
  integer :: maxmshd = 30

  logical :: linterface = .false.

  ! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
  ! *  NSPIND = 1                                                       *
  ! *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     *
  ! *          function, set up in the spin-independent non-relativstic *
  ! *          (l,m_l)-representation                                   *
  ! *                                                                   *
  ! *********************************************************************


end module global_variables


#ifdef CPP_MPI
!< MPI broadcast routine for global variables (i.e. array dimensions etc.)
subroutine Bcast_global_variables()
    use mpi
    use global_variables
    use mod_mympi,   only: master
    implicit none

    !< number of paramters that are broadcasted
    integer :: N

    !< blocklength of variuables in derived data type and list of MPI datatypes
    integer, allocatable  :: blocklen1(:), etype1(:)
    !< error status
    integer :: ierr
    !< derived data type for collective communication
    integer :: myMPItype1
    !< MPI addresses
    integer(kind=MPI_ADDRESS_KIND), allocatable :: disp1(:)
    !< base address of first entry
    integer(kind=MPI_ADDRESS_KIND) :: base

    N  = 59
    allocate(blocklen1(N), etype1(N), disp1(N), stat=ierr)
    if ( ierr/=0 ) stop 'error allocating arrays in bcast_global_variables'

    call MPI_Get_address(N,             disp1(1), ierr)
    call MPI_Get_address(irid,          disp1(2), ierr)
    call MPI_Get_address(krel,          disp1(3), ierr)
    call MPI_Get_address(nfund,         disp1(4), ierr)
    call MPI_Get_address(ipand,         disp1(5), ierr)   
    call MPI_Get_address(ngshd,         disp1(6), ierr)  
    call MPI_Get_address(ncleb,         disp1(7), ierr) 
    call MPI_Get_address(knoco,         disp1(8), ierr)
    call MPI_Get_address(iemxd,         disp1(9), ierr)
    call MPI_Get_address(irnsd,        disp1(10), ierr)
    call MPI_Get_address(nmaxd,        disp1(11), ierr)
    call MPI_Get_address(ishld,        disp1(12), ierr)
    call MPI_Get_address(naclsd,       disp1(13), ierr)   
    call MPI_Get_address(nspotd,       disp1(14), ierr)  
    call MPI_Get_address(ntperd,       disp1(15), ierr) 
    call MPI_Get_address(ntrefd,       disp1(16), ierr)
    call MPI_Get_address(nsheld,       disp1(17), ierr)
    call MPI_Get_address(ncelld,       disp1(18), ierr)
    call MPI_Get_address(nspind,       disp1(19), ierr)
    call MPI_Get_address(knosph,       disp1(20), ierr)
    call MPI_Get_address(korbit,       disp1(21), ierr)
    call MPI_Get_address(kpoibz,       disp1(22), ierr)
    call MPI_Get_address(wlength,      disp1(23), ierr) 
    call MPI_Get_address(nprincd,      disp1(24), ierr)
    call MPI_Get_address(nlayerd,      disp1(25), ierr) 
    call MPI_Get_address(natomimpd,    disp1(26), ierr)
    call MPI_Get_address(lmaxd,        disp1(27), ierr)
    call MPI_Get_address(lmmaxd,       disp1(28), ierr)
    call MPI_Get_address(lmgf0d,       disp1(29), ierr)
    call MPI_Get_address(alm,          disp1(30), ierr)
    call MPI_Get_address(almgf0,       disp1(31), ierr)
    call MPI_Get_address(ndim_slabinv, disp1(32), ierr)
    call MPI_Get_address(nembd,        disp1(33), ierr)
    call MPI_Get_address(nembd1,       disp1(34), ierr)
    call MPI_Get_address(nembd2,       disp1(35), ierr)
    call MPI_Get_address(nrd,          disp1(36), ierr)
    call MPI_Get_address(lm2d,         disp1(37), ierr)
    call MPI_Get_address(nclsd,        disp1(38), ierr)
    call MPI_Get_address(mmaxd,        disp1(39), ierr)
    call MPI_Get_address(npotd,        disp1(40), ierr)
    call MPI_Get_address(lmxspd,       disp1(41), ierr)
    call MPI_Get_address(lassld,       disp1(42), ierr)
    call MPI_Get_address(irmind,       disp1(43), ierr)
    call MPI_Get_address(nofgij,       disp1(44), ierr)
    call MPI_Get_address(nspindd,      disp1(45), ierr)
    call MPI_Get_address(nsatypd,      disp1(46), ierr)
    call MPI_Get_address(nrefd,        disp1(47), ierr)
    call MPI_Get_address(irmd,         disp1(48), ierr)
    call MPI_Get_address(naezd,        disp1(49), ierr)
    call MPI_Get_address(natypd,       disp1(50), ierr)
    call MPI_Get_address(lmpotd,       disp1(51), ierr)
    call MPI_Get_address(ntotd,        disp1(52), ierr)
    call MPI_Get_address(nrmaxd,       disp1(53), ierr)
    call MPI_Get_address(lmmaxso,      disp1(54), ierr)
    call MPI_Get_address(lpotd,        disp1(55), ierr)
    call MPI_Get_address(nchebd,       disp1(56), ierr)
    call MPI_Get_address(maxmshd,      disp1(57), ierr)
    call MPI_Get_address(linterface,   disp1(58), ierr)
    call MPI_Get_address(lnc,          disp1(59), ierr)

    ! find displacements of variables
    base  = disp1(1)
    disp1 = disp1 - base

    ! set length of variables in derived data type
    blocklen1(1:N) = 1

    ! set datatype of variables
    etype1(1:N-2) = MPI_INTEGER
    etype1(N-1:N) = MPI_LOGICAL

    ! create new Type structure for derived data type
    call MPI_Type_create_struct(N, blocklen1, disp1, etype1, myMPItype1, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in create_mpimask_t_inc'

    ! commit new type
    call MPI_Type_commit(myMPItype1, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error commiting create_mpimask_t_inc'

    ! broadcast derived data type
    
    call MPI_Bcast(N, 1, myMPItype1, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting t_inc'

    ! finally free auxiliary type and deallocate working arrays
    call MPI_Type_free(myMPItype1, ierr)
    deallocate(blocklen1, etype1, disp1, stat=ierr)
    if ( ierr/=0 ) stop 'error deallocating arrays in bcast_global_variables'

end subroutine Bcast_global_variables
#endif

