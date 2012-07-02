module KKRnano_mpi_mod
  implicit none

  SAVE

  integer :: my_SE_communicator   !< MPI communicator for (spin, energy)-group of process
  integer :: my_SE_comm_size      !< number of ranks in 'my_SE_communicator'
  integer :: my_SE_rank           !< rank of process in 'my_SE_communicator'
  logical :: is_Masterrank        !< true if process is the MASTERRANK
  integer :: my_spin_rank         !< The spin-index the process is working on
  integer :: my_energy_rank       !< The energy-group a process belongs to

  integer, private, parameter :: LMPID = 1 ! TODO: remove

  integer :: MYRANK
  integer :: LMPIC

  ! E-MPI
  integer :: EMPIC
  integer, private :: EMPIB

  ! S-MPI
  integer ::SMPIC
  integer, private ::SMPIB

  !     .. ACTV-MPI
  integer :: MYACTVRANK
  integer :: ACTVCOMM
  integer, private::ACTVGROUP

  ! ----------- arrays ----------------------------------------------------------

  ! -------------- MPI --------------------------------------------------------
  integer, dimension(:), allocatable, private :: MYLRANK
  integer, dimension(:), allocatable, private :: LCOMM
  integer, dimension(:), allocatable, private :: LGROUP
  integer, dimension(:), allocatable, private :: LSIZE

  !     .. S-MPI
  integer, dimension(:,:), allocatable, private :: SRANK
  integer, dimension(:,:), allocatable :: SMYRANK

  !     .. E-MPI
  integer, dimension(:,:), allocatable :: EMYRANK

  private :: fatalMemoryError
  private :: IMPI
  private :: allocate_KKRnano_mpi_arrays
  private :: deallocate_KKRnano_mpi_arrays

contains

  !> Display memory error message, stop program.
  !> @param msg An optional error message for display
  subroutine fatalMemoryError(msg)
    implicit none
    character(len=*), intent(in), optional :: msg

    write(*,*) "FATAL error, failure to (de)allocate memory."
    if (present(msg)) then
      write(*,*) msg
    end if

    stop

  end subroutine

  !------------------------------------------------------------------------------
  subroutine allocate_KKRnano_mpi_arrays(SMPID, EMPID, NAEZ)
    implicit none

    integer :: SMPID
    integer :: EMPID
    integer :: NAEZ

    integer :: memory_stat

    allocate(MYLRANK(LMPID*SMPID*EMPID), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("KKRnano_mpi_mod")
    MYLRANK = 0
    allocate(LCOMM(LMPID*SMPID*EMPID), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("KKRnano_mpi_mod")
    LCOMM = 0
    allocate(LGROUP(LMPID*SMPID*EMPID), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("KKRnano_mpi_mod")
    LGROUP = 0
    allocate(LSIZE(LMPID*SMPID*EMPID), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("KKRnano_mpi_mod")
    LSIZE = 0
    allocate(SRANK(SMPID,NAEZ*LMPID*EMPID), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("KKRnano_mpi_mod")
    SRANK = 0
    allocate(SMYRANK(SMPID,NAEZ*LMPID*EMPID), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("KKRnano_mpi_mod")
    SMYRANK = 0
    allocate(EMYRANK(EMPID,NAEZ*LMPID*SMPID), stat = memory_stat)
    if(memory_stat /= 0) call fatalMemoryError("KKRnano_mpi_mod")
    EMYRANK = 0
  end subroutine

  !------------------------------------------------------------------------------
  subroutine deallocate_KKRnano_mpi_arrays()
    implicit none

    integer :: memory_stat

    deallocate(MYLRANK, stat = memory_stat)
    deallocate(LCOMM, stat = memory_stat)
    deallocate(LGROUP, stat = memory_stat)
    deallocate(LSIZE, stat = memory_stat)
    deallocate(SRANK, stat = memory_stat)
    deallocate(SMYRANK, stat = memory_stat)
    deallocate(EMYRANK, stat = memory_stat)
  end subroutine

  !------------------------------------------------------------------------------
  subroutine initialiseKKRnano_mpi_com(SMPID, EMPID, NAEZ, nthrds)
    implicit none

    integer :: SMPID
    integer :: EMPID
    integer :: NAEZ
    integer :: nthrds

    !------- local

      ! not needed outside of this routine
    integer :: NROFNODES
    integer :: ACTVSIZE

    call allocate_KKRnano_mpi_arrays(SMPID, EMPID, NAEZ)

    call IMPI(NAEZ,MYRANK,NROFNODES, &
    LMPIC,MYLRANK,LGROUP,LCOMM,LSIZE, &
    SMPIB,SMPIC,SRANK,SMYRANK, &
    EMPIB,EMPIC,EMYRANK, &
    MYACTVRANK,ACTVGROUP,ACTVCOMM,ACTVSIZE, &
    lmpid, smpid, empid, nthrds)

    ! Here comes the initialisation of the variables that define the status of
    ! a process
    if (LMPIC > 0) then
      my_SE_communicator = LCOMM(LMPIC)
      my_SE_comm_size = LSIZE(LMPIC)
      my_SE_rank = MYLRANK(LMPIC)
      my_spin_rank = SRANK(SMPIB,SMPIC)
      my_energy_rank = EMPIB
      is_Masterrank = (MYLRANK(1) == 0)
    else
      ! process inactive
      my_SE_communicator = -1
      my_SE_comm_size = -1
      my_SE_rank = -1
      my_spin_rank = -1
      my_energy_rank = -1
      is_Masterrank = .false.
    end if

  end subroutine

  !------------------------------------------------------------------------------
  subroutine finaliseKKRnano_mpi_com()
    implicit none

    include 'mpif.h'

    integer :: IERR

    ! Free communicators and groups ..

    ! only active processes
    if (LMPIC /= 0) then !if (MYACTVRANK >= 0) then

      !if (MYLRANK(LMPIC)>=0) then
      call MPI_COMM_FREE(LCOMM(LMPIC),IERR)
      !endif

      call MPI_GROUP_FREE(LGROUP(LMPIC),IERR)

      call MPI_COMM_FREE(ACTVCOMM,IERR)
      call MPI_GROUP_FREE(ACTVGROUP,IERR)

    end if

    ! all processes
    call deallocate_KKRnano_mpi_arrays()

    call MPI_BARRIER(MPI_COMM_WORLD,IERR)
    call MPI_FINALIZE(IERR)

  end subroutine

  !------------------------------------------------------------------------------
  subroutine IMPI( &
  NAEZ, &                                     ! > in
  MYRANK,NROFNODES, &                         ! < out
  LMPIC,MYLRANK,LGROUP,LCOMM,LSIZE, &         ! < out
  SMPIB,SMPIC,SRANK,SMYRANK, &                ! < out
  EMPIB,EMPIC,EMYRANK, &                ! < out
  MYACTVRANK,ACTVGROUP,ACTVCOMM,ACTVSIZE, &   ! < out
  lmpid, smpid, empid, nthrds)                ! < in  ! new input parameters after inc.p removal

    ! ======================================================================
    !                       build MPI_Groups
    ! ======================================================================

    ! the required information for l-parallelization are
    ! intitialized.

    ! all depends on the additional switch in inc.p:

    ! MPI_L = 1  .... still one process per atom
    ! MPI_L > 1  .... LMPI number of processes running per atom

    ! parallelization strategy:
    ! only the solution of the Dyson equation, resp. the inversion of
    ! the Lloyd-matrix are parallelized of LM. All tasks before are
    ! executed anyway, to save communication. After INVIT and LLYINVIT
    ! finished, the results for different L are communicated to the
    ! root process by MPI_REDUCE(...,MPI_SUM,...) and from here on only
    ! the rout processes are running.
    !                               Alexander Thiess, 19th of November 2009

    !                               f90 conversion: Elias Rabel, Oct 2011
    ! =======================================================================

    !use mpi
    implicit none
  include 'mpif.h'

    integer, intent(in) :: lmpid
    integer, intent(in) :: smpid
    integer, intent(in) :: empid
    integer, intent(in) :: nthrds

    integer::I1
    integer::NAEZ
    integer::ISPIN
    integer::IRANK
    !     .. MPI ..
    integer::WGROUP
    integer::IERR
    !     .. N-MPI
    integer::MYRANK
    integer::NROFNODES
    !     .. L-MPI

    integer::MYLRANK    (LMPID*SMPID*EMPID)
    integer::LRANKS(NAEZ,LMPID*SMPID*EMPID)
    integer::LCOMM      (LMPID*SMPID*EMPID)
    integer::LGROUP     (LMPID*SMPID*EMPID)
    integer::LSIZE      (LMPID*SMPID*EMPID)
    integer::LMPI
    integer::LMPIC

    !     .. S-MPI
    integer::SMYRANK(SMPID,NAEZ*LMPID*EMPID)
    integer::SRANK  (SMPID,NAEZ*LMPID*EMPID)
    integer::SMPIC
    integer::SMPIB

    !     .. E-MPI
    integer::EMYRANK(EMPID,NAEZ*LMPID*SMPID)
    integer::EMPI
    integer::EMPIC
    integer::EMPIB

    !     .. ACTV-MPI
    integer::MYACTVRANK
    integer::ACTVRANKS(NAEZ*LMPID*SMPID*EMPID)
    integer::ACTVCOMM
    integer::ACTVGROUP
    integer::ACTVSIZE

    !$ integer::omp_get_max_threads

    ! intialize MPI and standard groups and communicator

    call MPI_INIT(IERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,NROFNODES,IERR)

    ! define ranks for all non-standard groups of processors

    LMPIC = 0

    SMPIC = 0
    SMPIB = 0

    EMPIC = 0
    EMPIB = 0


    do I1 = 1, NAEZ
      do ISPIN = 1, SMPID
        do LMPI = 1, LMPID
          do EMPI = 1, EMPID

            IRANK = (ISPIN-1)*(NAEZ*LMPID*EMPID) + &
            (LMPI-1)*NAEZ*EMPID + &
            (I1-1)*EMPID + &
            EMPI - 1

            LRANKS(I1,(ISPIN-1)*LMPID*EMPID+(LMPI-1)*EMPID+EMPI) = IRANK

            SMYRANK(ISPIN,(LMPI-1)*NAEZ*EMPID+(I1-1)*EMPID+EMPI) = IRANK
            SRANK(ISPIN,(LMPI-1)*NAEZ*EMPID+(I1-1)*EMPID+EMPI)   = ISPIN-1

            EMYRANK(EMPI,(ISPIN-1)*NAEZ*LMPID+(I1-1)*LMPID+LMPI)= IRANK
            !ERANK(EMPI,(ISPIN-1)*NAEZ*LMPID+(I1-1)*LMPID+LMPI)  = LMPI-1

            ACTVRANKS(IRANK+1) = IRANK

            if (MYRANK == IRANK) then

              LMPIC  = (ISPIN-1)*LMPID*EMPID+(LMPI-1)*EMPID + EMPI

              SMPIC  = (LMPI-1)*NAEZ*EMPID  +(I1-1)*EMPID   + EMPI
              SMPIB  = ISPIN

              EMPIC  = (ISPIN-1)*NAEZ*LMPID +(I1-1)*LMPID   + LMPI
              EMPIB  = EMPI

            endif

          enddo
        enddo
      enddo
    enddo

    ! build groups and communicators


    ! ACTIVE GROUP (ACTVGROUP) ...............................................

    call MPI_COMM_GROUP(MPI_COMM_WORLD,WGROUP,IERR) ! get group for MPI_COMM_WORLD
    call MPI_GROUP_INCL(WGROUP,NAEZ*LMPID*SMPID*EMPID, &
    ACTVRANKS(1), &
    ACTVGROUP,IERR) !create a group ACTVGROUP
    call MPI_COMM_CREATE(MPI_COMM_WORLD,ACTVGROUP,ACTVCOMM,IERR)

    MYACTVRANK=-2
    call MPI_GROUP_RANK(ACTVGROUP,MYACTVRANK,IERR)
    call MPI_GROUP_SIZE(ACTVGROUP,ACTVSIZE,IERR)

    ! LGROUP ...............................................................

    do LMPI=1,LMPID*SMPID*EMPID
      call MPI_COMM_GROUP(MPI_COMM_WORLD, WGROUP, IERR) ! get default
      call MPI_GROUP_INCL(WGROUP, NAEZ, LRANKS(1,LMPI), &
      LGROUP(LMPI), IERR) !create groups LGROUP
      call MPI_COMM_CREATE(MPI_COMM_WORLD, LGROUP(LMPI), &
      LCOMM(LMPI), IERR) !create communicators for
    enddo

    do LMPI=1,LMPID*SMPID*EMPID
      MYLRANK(LMPI) = -2
      call MPI_GROUP_RANK(LGROUP(LMPI),MYLRANK(LMPI),IERR)
      call MPI_GROUP_SIZE(LGROUP(LMPI),LSIZE(LMPI),IERR)
    enddo

    !.......................................................................

    if (MYRANK == 0) then
      write(6,'(79(1H=))')
      write(6,*) '  total no. of MPI ranks  = ',NROFNODES
      !$omp parallel
      !$omp single
      !$  write(6,*) '  OMP max. threads        = ',omp_get_max_threads()
      !$omp end single
      !$omp end parallel
      write(6,*) '  groups of processes created'
      write(6,*) '  NMPI                    = ',NAEZ
      write(6,*) '  LMPI                    = ',LMPID
      write(6,*) '  SMPI                    = ',SMPID
      write(6,*) '  EMPI                    = ',EMPID
      write(6,*) '  NTHRDS                  = ',NTHRDS
      write(6,*) '  MPI-processes (active)  = ',NAEZ*LMPID*SMPID*EMPID
      write(6,*) '  total no. of tasks      = ',NAEZ*LMPID*SMPID*EMPID*NTHRDS
      write(6,'(79(1H=))')
    endif


  end subroutine IMPI

end module KKRnano_mpi_mod
