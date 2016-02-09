#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

#define COMMCHECK(IERR) if((IERR) /= 0) then; write(*,*) "ERROR: Communication failure: ", __FILE__, __LINE__; STOP; endif


module EBalanceHandler_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: EBalanceHandler, destroy !, destroyEBalanceHandler
  public :: createEBalanceHandler, startEBalanceTiming, stopEBalanceTiming
  public :: initEBalanceHandler, updateEBalance_com, setEqualDistribution
  
  type EBalanceHandler
    integer, allocatable :: eproc(:)
    integer, allocatable :: eproc_old(:)
    real, allocatable :: etime(:)
    integer :: ierlast
    integer :: num_eprocs_empid
    real :: TIME_E
    real :: TIME_EX
    logical :: equal_distribution
  endtype
  
!   interface create
!     module procedure createEBalanceHandler
!   endinterface
!
!   interface init
!     module procedure initEBalanceHandler
!   endinterface
 
  interface destroy
    module procedure destroyEBalanceHandler
  endinterface

  contains

  !----------------------------------------------------------------------------
  !> Creates EBalanceHandler.
  !> Call initEBalanceHandler before use
  subroutine createEBalanceHandler(balance, num_epoints)
    type(EBalanceHandler), intent(inout) :: balance
    integer, intent(in) :: num_epoints

    integer :: memory_stat

    ALLOCATECHECK(balance%eproc(num_epoints))
    ALLOCATECHECK(balance%eproc_old(num_epoints))
    ALLOCATECHECK(balance%etime(num_epoints))

    balance%eproc = 0
    balance%eproc_old = 0
    balance%etime = 0.0d0

    balance%ierlast = num_epoints
    balance%num_eprocs_empid = 0
    balance%equal_distribution = .false.

  endsubroutine ! create

  !----------------------------------------------------------------------------
  !> Initialises EBalanceHandler.
  subroutine initEBalanceHandler(balance, mp)
    use KKRnanoParallel_mod, only: KKRnanoParallel
    type(EBalanceHandler), intent(inout) :: balance
    type(KKRnanoParallel), intent(in) :: mp

    integer :: file_points, file_procs, ios
    integer, parameter :: FILEHANDLE = 50
    character(len=*), parameter :: filename = 'ebalance'

    balance%num_eprocs_empid = mp%numEnergyRanks

    if (balance%num_eprocs_empid == 1) then
      balance%eproc = 1
      balance%eproc_old = 1
      return
    endif

    if (mp%myWorldRank == 0) then
      ! TODO: check and read ebalance file

      call EBALANCE1(balance%ierlast, balance%eproc, balance%eproc_old, &
                     balance%num_eprocs_empid, balance%ierlast)

      open(FILEHANDLE, file=filename, form='formatted', status='old', action='read', iostat=ios)
      if (ios == 0) then
        call readEBalanceHeader(filehandle, file_points, file_procs)

        if (balance%ierlast == file_points .and. balance%num_eprocs_empid == file_procs) then

          call readEBalanceDistribution(filehandle, balance%eproc, file_points)

        else
          warn(6, "Bad file "-filename-" provided")
          write(*,*) "WARNING: Bad ebalance file provided."
        endif
        close(FILEHANDLE)
      endif

    endif

    ! bcast ebalance distribution
    call bcastEBalance_com(balance, mp)

    balance%eproc_old = balance%eproc

  endsubroutine ! init

  !----------------------------------------------------------------------------
  !> Sets option for equal work distribution.
  !> Equal work distribution is used for DOS calculation.
  !> Default is: no equal work distribution
  subroutine setEqualDistribution(balance, flag)
    type(EBalanceHandler), intent(inout) :: balance
    logical, intent(in) :: flag

    balance%equal_distribution = flag

  endsubroutine ! set

  !----------------------------------------------------------------------------
  !> (Re)Starts measurement for dynamical load balancing
  subroutine startEBalanceTiming(balance, ie)
    type(EBalanceHandler), intent(inout) :: balance
    integer, intent(in) :: ie

    ! Start timing:
    ! save initial time of calculation
    call CPU_TIME(balance%time_e)
    balance%etime(ie) = 0.0d0

  endsubroutine

  !----------------------------------------------------------------------------
  !> Finishes time measurement for energy-point ie
  subroutine stopEBalanceTiming(balance, ie)
    type(EBalanceHandler), intent(inout) :: balance
    integer, intent(in) :: ie

    real :: time_finish

    call CPU_TIME(time_finish)

    balance%etime(ie) = time_finish - balance%time_e

  endsubroutine

  !----------------------------------------------------------------------------
  !> Gathers all timings and redistributes energy processes.
  subroutine updateEBalance_com(balance, mp)
    use KKRnanoParallel_mod, only: KKRnanoParallel
    implicit none
    include 'mpif.h'

    type(EBalanceHandler), intent(inout) :: balance
    type(KKRnanoParallel), intent(in) :: mp

    integer :: npnt1, ierr

    real MTIME(balance%ierlast)

    npnt1 = 1; if (balance%equal_distribution) npnt1 = 0

    ! TODO: extract allreduce -> reduce
    ! TODO: let only master rank calculate
    ! TODO: broadcast
    if (balance%num_eprocs_empid > 1) then

      call MPI_REDUCE(balance%ETIME,MTIME,balance%ierlast,MPI_REAL,MPI_MAX, 0, mp%myActiveComm,IERR)

      COMMCHECK(IERR)

      ! save old rank distribution
      balance%eproc_old = balance%eproc

      ! only Masterrank calculates new work distribution
      if (mp%myWorldRank == 0) then
        call EBALANCE2(balance%ierlast,npnt1, mp%myWorldRank, mp%myActiveComm, MTIME,balance%EPROC,balance%EPROC_old, balance%num_eprocs_empid, balance%ierlast)
      endif

      call bcastEBalance_com(balance, mp)

   endif

  endsubroutine

  !----------------------------------------------------------------------------
  !> Destroys EBalanceHandler.
  !> Call initEBalanceHandler before use
  subroutine destroyEBalanceHandler(balance)
    type(EBalanceHandler), intent(inout) :: balance

    integer :: memory_stat

    DEALLOCATECHECK(balance%eproc)
    DEALLOCATECHECK(balance%eproc_old)
    DEALLOCATECHECK(balance%etime)
  endsubroutine


!==============================================================================
! H E L P E R       R O U T I N E S
!==============================================================================

  !----------------------------------------------------------------------------
  !> Broadcast ebalance information from rank 0 to all ranks.
  subroutine bcastEBalance_com(balance, mp)
    use KKRnanoParallel_mod, only: KKRnanoParallel
    include 'mpif.h'

    type(EBalanceHandler), intent(inout) :: balance
    type(KKRnanoParallel), intent(in) :: mp

    integer :: ierr

    ! save old ebalance information WRONG!!!! already done
    ! balance%eproc_old = balance%eproc

    call MPI_BCAST(balance%eproc, balance%ierlast, MPI_INTEGER, 0, mp%myActiveComm, ierr)

    COMMCHECK(IERR)

  endsubroutine ! bcast

  !---------------------------------------------------------------------------
  !> Reads header of opened ebalance file and returns number of points
  !> and process-groups from file
  subroutine readEBalanceHeader(filehandle, npoints, nprocs)
    integer, intent(in) :: filehandle
    integer, intent(out) :: npoints
    integer, intent(out) :: nprocs

    !skip first 3 lines - those are comments
    read(filehandle, *)
    read(filehandle, *)
    read(filehandle, *)
    read(filehandle, *) npoints, nprocs
  endsubroutine ! read

  !---------------------------------------------------------------------------
  !> Reads body of opened ebalance file and returns energy group distribution
  !> note: no proper error handling
  subroutine readEBalanceDistribution(filehandle, eproc, num_points)
    integer, intent(in) :: filehandle
    integer, intent(inout) :: eproc(:)
    integer, intent(in) :: num_points

    integer :: idummy, ie
    real :: rdummy

    do ie = 1, num_points
      read(filehandle, *) idummy, eproc(ie), rdummy
    enddo ! ie
    
  endsubroutine ! read

!------------------------------------------------------------------------------
!                                                optional: equal
subroutine ebalance1(ierlast, eproc, eproco, empid, iemxd, equal)
  ! ======================================================================
  !                       build mpi_groups
  ! ======================================================================

  ! this routine balances the distribution of energy points
  ! to the set of processors of empi=1,...,empid


  !                               alexander thiess, 7th of december 2009
  ! =======================================================================

  include 'mpif.h'

  integer, intent(in) :: ierlast, empid, iemxd ! iemxd is to be removed
  integer, intent(out) :: eproc(iemxd), eproco(iemxd)
  logical, intent(in), optional :: equal

  integer :: ier,ieri
  logical :: flag

  if (present(equal)) then
     flag = equal
  else
     flag = .false.
  endif

  !=======================================================================
  ! 1st guess >>>
  !=======================================================================
  !     trying to figure out what this code means:
  !     *) if there is only one energy process -> it calculates everything
  !     *) if there are two processes -> last two points calculated by
  !        process 1, rest is calculated by process 2
  !     *) 3 processes: 1st process calculates only last point
  !                     2nd process calculates first 2/3 of all points
  !                     3rd process calculates last 1/3 of points except
  !                     last
  !     *) more than 3 processes: 1st proc. last point
  !                               2nd proc. first 1/2 of all points
  !                               other procs.: rest


  if (empid == 1) then
    do ier = 1, ierlast
      eproc(ier) = 1
    enddo ! ier
  elseif (empid == 2) then
    eproc(ierlast)   = 1
    eproc(ierlast-1) = 1
    do ier = 1, ierlast-2
      eproc(ier)     = 2
    enddo ! ier
  elseif (empid == 3) then
    eproc(ierlast)   = 1
    ieri = floor(ierlast/3.*2.)
    do ier = 1, ieri
      eproc(ier) = 2
    enddo
    do ier = ieri+1, ierlast-1
      eproc(ier) = 3
    enddo
  else
    eproc(ierlast) = 1
    ieri = ierlast / 2
    do ier = 1, ieri
      eproc(ier) = 2
    enddo ! ier
    do ier = ieri+1, ierlast-1
      eproc(ier) = mod(ier, empid-2) + 3
    enddo ! ier
  endif

  !     if dos shall be calculated make use of special scheme allowing for
  !     broader energy-parallelization

  ! remark e.r.:
  ! the amazing "special"-scheme just means equal work distribution

  !if (test('dos     ')) then
  if (flag) then
    do ier = 1, ierlast
      eproc(ier) = mod(ier,empid) + 1
    enddo ! ier
  endif

  eproco = eproc

  !=======================================================================
  ! >>> 1st guess
  !=======================================================================
endsubroutine ! ebalance1

!==============================================================================

subroutine ebalance2(ierlast, npnt1, myactvrank, activecomm, mtime, eproc, eproco, empid, iemxd)
  ! ======================================================================
  !                       build mpi_groups
  ! ======================================================================
  ! this routine balances the distribution of energy points
  ! to the set of processors of empi=1,...,empid
  !                               alexander thiess, 7th of december 2009
  ! =======================================================================
  integer empid
  integer iemxd
  integer        ierlast,npnt1
  integer        eproc(iemxd), eproco(iemxd), epoint(iemxd,empid)
  real           aimtm,sumtm
  integer empi,ier,ieri

  real           proctm(empid), mtime(iemxd)
  integer        myactvrank, activecomm

  !=======================================================================
  ! use timing of iter-1 >>>
  !=======================================================================
  !     gather on all processors all timings
  !call mpi_allreduce(etime,mtime,iemxd,mpi_real,mpi_max, activecomm,ierr)


  !----------------------------------------------------------------------
  !     analyze balance of last iteration
  !----------------------------------------------------------------------

  do empi = 1, empid
    proctm(empi) = 0.d0
  enddo ! empi
  
  sumtm = 0.d0
  do ier = 1, ierlast
    sumtm = sumtm + mtime(ier)
    proctm(eproc(ier))  = proctm(eproc(ier)) + mtime(ier)
    eproco(ier) = eproc(ier)
  enddo ! ier

  do empi = 1, empid
    if (myactvrank == 0) write(6,*) 'EMPID: group ',empi,' took',proctm(empi),'%:',100.*proctm(empi)/sumtm
  enddo ! empi

  !----------------------------------------------------------------------
  !     intitialize
  !----------------------------------------------------------------------

  sumtm = 0.0d0

  do ier = 1, ierlast
    sumtm = sumtm + mtime(ier)
    do empi = 1, empid
      epoint(ier,empi) = 0
    enddo ! empi
    eproc(ier) = 0
  enddo ! ier

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------


  !     goal for total time on each processor 1,...,empid

  aimtm = sumtm/empid

  !     1st processor (empi=1) acts always on ier=ierlast
  !     and if more time is availible on energy points starting with npnt1
  !     2nd ... empid processors sum up timing from ier=1 up to ierlast-1
  !     taking time limits into account

  !     loop over processors empi=1,empid >>>
  do empi = 1, empid

    proctm(empi) = 0.d0

    if (empi == 1) then
      proctm(1) = mtime(ierlast)
      epoint(ierlast,1) = ierlast
      eproc(ierlast) = 1
      ieri = npnt1
    else
      ieri = 1
    endif

    if (empid == 1) ieri = 1


    !     quick fix: if npnt1=0 then ieri=0 -> leads to read out of bounds
    !     e.r.
    if (ieri == 0) ieri = 1

    !     loop over processors ier=ieri,ierlast >>>
    do ier = ieri, ierlast
      if (empi < empid) then
        if ((proctm(empi) < aimtm) .and. (eproc(ier) == 0)) then
          proctm(empi) = proctm(empi) + mtime(ier)
          epoint(ier,empi)= ier
          eproc(ier) = empi
        endif
      elseif (empi == empid) then
        if (eproc(ier) == 0) then
          proctm(empi) = proctm(empi) + mtime(ier)
          epoint(ier,empi)= ier
          eproc(ier) = empi
        endif
      endif
    enddo
    !     >>> loop over processors ier=ieri,ierlast

    if (myactvrank == 0) write(6,*) 'EMPID: group ',empi,'will take',proctm(empi),'%:',100.*proctm(empi)/sumtm

  enddo ! empi
  !     >>> loop over processors empi=1,empid


!   write information on load-balancing to formatted file 'balance'

!   inquire(file='stop', exist=stopit)
!   if (iter == scfsteps) then
    if (myactvrank == 0) then
      open (50, file='ebalance', form='formatted', status='replace', action='write')
      write(50, *) "# Energy load-balancing file"
      write(50, *) "# 1st line: number of E-points, number of E-processes."
      write(50, *) "# point --- process --- timing used"
      write(50, *) ierlast,empid

      do ier = 1, ierlast
        write(50,*) ier, eproc(ier), mtime(ier)
      enddo ! ier
      
      close(50)
    endif
!   endif

endsubroutine ! ebalance2


endmodule ! EBalanceHandler_mod
