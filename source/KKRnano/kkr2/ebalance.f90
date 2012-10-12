module EBalanceHandler_mod
  implicit none

  type EBalanceHandler
    PRIVATE
    integer, dimension(:), allocatable :: eproc
    integer, dimension(:), allocatable :: eproc_old
    real, dimension(:), allocatable :: etime
    integer :: ierlast
    integer :: num_eprocs_empid
    real::TIME_E
    real::TIME_EX
  end type

  CONTAINS

  !----------------------------------------------------------------------------
  !> Creates EBalanceHandler.
  !> Call initEBalanceHandler before use
  subroutine createEBalanceHandler(balance, num_epoints)
    implicit none
    type (EBalanceHandler), intent(inout) :: balance
    integer, intent(in) :: num_epoints

    allocate(balance%eproc(num_epoints))
    allocate(balance%eproc_old(num_epoints))
    allocate(balance%etime(num_epoints))

    balance%eproc = 0
    balance%eproc_old = 0
    balance%etime = 0.0d0

    balance%ierlast = num_epoints
    balance%num_eprocs_empid = 0

  end subroutine

  !----------------------------------------------------------------------------
  !> Initialises EBalanceHandler.
  subroutine initEBalanceHandler(balance, my_mpi)
    use KKRnanoParallel_mod
    implicit none
    type (EBalanceHandler), intent(inout) :: balance
    type (KKRnanoParallel), intent(in) :: my_mpi

    balance%num_eprocs_empid = getNumEnergyRanks(my_mpi)

    if (isMasterRank(my_mpi)) then
      ! TODO: check and read ebalance file

      call EBALANCE1(balance%ierlast, balance%eproc, balance%eproc_old, &
                     balance%num_eprocs_empid, balance%ierlast)

    end if

  end subroutine

  !----------------------------------------------------------------------------
  !> Broadcast ebalance information to all ranks.
  subroutine bcastEBalance_com(balance, my_mpi)
    use KKRnanoParallel_mod
    implicit none

    include 'mpif.h'

    type (EBalanceHandler), intent(inout) :: balance
    type (KKRnanoParallel), intent(in) :: my_mpi

    !-----
    integer :: ierr

    ! save old ebalance information
    balance%eproc_old = balance%eproc

    call MPI_BCAST(balance%eproc,balance%ierlast,MPI_INTEGER, &
    getMasterRank(my_mpi),getMyActiveCommunicator(my_mpi), ierr)

  end subroutine

  !----------------------------------------------------------------------------
  !> (Re)Starts measurement for dynamical load balancing
  subroutine startEBalanceTiming(balance, ie)
    use KKRnanoParallel_mod

    implicit none

    type (EBalanceHandler), intent(inout) :: balance
    integer, intent(in) :: ie

    ! Start timing:
    ! save initial time of calculation
    call CPU_TIME(balance%time_e)
    balance%etime(ie) = 0.0d0

  end subroutine

  !----------------------------------------------------------------------------
  !> Finishes time measurement for energy-point ie
  subroutine stopEBalanceTiming(balance, ie)
    use KKRnanoParallel_mod
    implicit none

    type (EBalanceHandler), intent(inout) :: balance
    integer, intent(in) :: ie

    real :: time_finish

    call CPU_TIME(time_finish)

    balance%etime(ie) = time_finish - balance%time_e

  end subroutine

  !----------------------------------------------------------------------------
  ! Gathers all timings and redistributes energy processes.
  subroutine updateEBalance_com(balance, my_mpi)
    use KKRnanoParallel_mod
    implicit none

    type (EBalanceHandler), intent(inout) :: balance
    type (KKRnanoParallel), intent(in) :: my_mpi

    ! TODO TODO TODO
    ! CAST ETIME to real!!!
    ! call EBALANCE2   GET NPNT1

    integer :: NPNT1

    NPNT1 = 1

    call EBALANCE2(balance%ierlast,NPNT1, getMyWorldRank(my_mpi), &
    getMyActiveCommunicator(my_mpi), &
    balance%ETIME,balance%EPROC,balance%EPROC_old, &
    balance%num_eprocs_empid, balance%ierlast)

  end subroutine

  !----------------------------------------------------------------------------
  !> Destroys EBalanceHandler.
  !> Call initEBalanceHandler before use
  subroutine destroyEBalanceHandler(balance)
    use KKRnanoParallel_mod
    implicit none
    type (EBalanceHandler), intent(inout) :: balance

    deallocate(balance%eproc)
    deallocate(balance%eproc_old)
    deallocate(balance%etime)

  end subroutine

end module

!------------------------------------------------------------------------------
subroutine EBALANCE1(IERLAST, EPROC, EPROCO, empid, iemxd)
  ! ======================================================================
  !                       build MPI_Groups
  ! ======================================================================

  ! this routine balances the distribution of energy points
  ! to the set of processors of EMPI=1,...,EMPID


  !                               Alexander Thiess, 7th of December 2009
  ! =======================================================================

  implicit none

  include 'mpif.h'

  integer empid
  integer iemxd

  !     .. global scalars ..
  integer IERLAST

  integer        EPROC(IEMXD), EPROCO(IEMXD)

  integer        IER,IERI
  logical        TEST

  !=======================================================================
  ! 1st guess >>>
  !=======================================================================
  !     Trying to figure out what this code means:
  !     *) if there is only one energy process -> it calculates everything
  !     *) if there are two processes -> last two points calculated by
  !        process 1, rest is calculated by process 2
  !     *) 3 processes: 1st process calculates only last point
  !                     2nd process calculates first 2/3 of all points
  !                     3rd process calculates last 1/3 of points except
  !                     last
  !     *) MORE than 3 processes: 1st proc. last point
  !                               2nd proc. first 1/2 of all points
  !                               other procs.: rest


  if (EMPID.eq.1) then
    do IER=1,IERLAST
      EPROC(IER)     = 1
    enddo
  elseif (EMPID.eq.2) then
    EPROC(IERLAST)   = 1
    EPROC(IERLAST-1) = 1
    do IER=1,(IERLAST-2)
      EPROC(IER)     = 2
    enddo
  elseif (EMPID.eq.3) then
    EPROC(IERLAST)   = 1
    IERI = FLOOR(REAL(IERLAST)/(REAL(3))*REAL(2))
    do IER=1, IERI
      EPROC(IER) = 2
    enddo
    do IER=(IERI+1), (IERLAST-1)
      EPROC(IER) = 3
    enddo
  else
    EPROC(IERLAST)   = 1
    IERI = IERLAST / 2
    do IER=1, IERI
      EPROC(IER) = 2
    enddo
    do IER=(IERI+1), (IERLAST-1)
      EPROC(IER) = mod(IER, (EMPID-2)) + 3
    enddo
  endif

  !     if DOS shall be calculated make use of special scheme allowing for
  !     broader Energy-parallelization

  if (TEST('DOS     ')) then
    do IER=1,IERLAST
      EPROC(IER)     = MOD(IER,EMPID) + 1
    enddo
  endif

  EPROCO = EPROC


  !     check whether information on energy-point load-balancing is
  !     available

  !     INQUIRE(FILE='balance',EXIST=READIT)

  !       IF ((ITER.EQ.1).AND.(READIT)) THEN

  !         OPEN (50,FILE='balance',FORM='unformatted')

  !         READ (50,IOSTAT=IOS) BLNCD,BLNCD1

  !         IF (IOS.EQ.0.AND.BLNCD.EQ.IERLAST.AND.BLNCD1.EQ.EMPID) THEN
  !           READ (50) EPROC
  !         ENDIF

  !         CLOSE(50)

  !       ENDIF

  !=======================================================================
  ! >>> 1st guess
  !=======================================================================
end subroutine

!==============================================================================

subroutine EBALANCE2(IERLAST,NPNT1, MYACTVRANK,ACTVCOMM, &
ETIME,EPROC,EPROCO, &
empid, iemxd)
  ! ======================================================================
  !                       build MPI_Groups
  ! ======================================================================

  ! this routine balances the distribution of energy points
  ! to the set of processors of EMPI=1,...,EMPID


  !                               Alexander Thiess, 7th of December 2009
  ! =======================================================================

  implicit none

  include 'mpif.h'

  integer empid
  integer iemxd

  integer        IERLAST,NPNT1

  real           ETIME(IEMXD)
  integer        EPROC(IEMXD), EPROCO(IEMXD), EPOINT(IEMXD,EMPID)

  !     .. local scalars ..
  real           AIMTM,SUMTM
  integer EMPI,IER,IERI


  !     .. local arrays ..
  real           PROCTM(EMPID), MTIME(IEMXD)

  !     .. MPI ACTV ..
  integer        MYACTVRANK,ACTVCOMM,IERR

  external       LSAME

  !=======================================================================
  ! use timing of ITER-1 >>>
  !=======================================================================

  !     gather on all processors all timings

  call MPI_ALLREDUCE(ETIME,MTIME,IEMXD,MPI_REAL,MPI_MAX, &
  ACTVCOMM,IERR)


  !----------------------------------------------------------------------
  !     analyze balance of last iteration
  !----------------------------------------------------------------------

  do EMPI=1,EMPID
    PROCTM(EMPI)        = 0.0D0
  enddo

  SUMTM                 = 0.0D0

  do IER=1,IERLAST
    SUMTM               = SUMTM + MTIME(IER)
    PROCTM(EPROC(IER))  = PROCTM(EPROC(IER)) + MTIME(IER)

    EPROCO(IER)         = EPROC(IER)
  enddo

  do EMPI=1,EMPID
    if (MYACTVRANK.eq.0) &
    write(6,*) 'EMPID: group ',EMPI,' took',PROCTM(EMPI),'%:', &
    REAL(100)*PROCTM(EMPI)/SUMTM
  enddo

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------


  !----------------------------------------------------------------------
  !     intitialize
  !----------------------------------------------------------------------

  SUMTM = 0.0D0

  do IER = 1, IERLAST
    SUMTM         = SUMTM + MTIME(IER)
    do EMPI = 1, EMPID
      EPOINT(IER,EMPI) = 0
    enddo
    EPROC(IER)    = 0
  enddo

  !----------------------------------------------------------------------
  !----------------------------------------------------------------------


  !     goal for total time on each processor 1,...,EMPID

  AIMTM = SUMTM/EMPID

  !     1st processor (EMPI=1) acts always on IER=IERLAST
  !     and if more time is availible on energy points starting with NPNT1
  !     2nd ... EMPID processors sum up timing from IER=1 up to IERLAST-1
  !     taking time limits into account

  !     loop over processors EMPI=1,EMPID >>>
  do EMPI=1,EMPID

    PROCTM(EMPI)        = 0.0D0

    if (EMPI.eq.1) then
      PROCTM(1)         = MTIME(IERLAST)
      EPOINT(IERLAST,1) = IERLAST
      EPROC(IERLAST)    = 1
      IERI              = NPNT1
    else
      IERI              = 1
    endif

    if (EMPID.eq.1) then
      IERI              = 1
    endif


    !     Quick fix: if NPNT1=0 then IERI=0 -> leads to read out of bounds
    !     E.R.
    if (IERI .eq. 0) then
      IERI = 1
    end if

    !     loop over processors IER=IERI,IERLAST >>>
    do IER = IERI, IERLAST
      if (EMPI.lt.EMPID) then
        if ((PROCTM(EMPI).lt.AIMTM).and.(EPROC(IER).eq.0)) then
          PROCTM(EMPI)    = PROCTM(EMPI) + MTIME(IER)
          EPOINT(IER,EMPI)= IER
          EPROC(IER)      = EMPI
        endif
      elseif (EMPI.eq.EMPID) then
        if (EPROC(IER).eq.0) then
          PROCTM(EMPI)    = PROCTM(EMPI) + MTIME(IER)
          EPOINT(IER,EMPI)= IER
          EPROC(IER)      = EMPI
        endif
      endif
    enddo
    !     >>> loop over processors IER=IERI,IERLAST

    if (MYACTVRANK.eq.0) &
    write(6,*) 'EMPID: group ',EMPI,'will take',PROCTM(EMPI),'%:', &
    REAL(100)*PROCTM(EMPI)/SUMTM


  enddo
  !     >>> loop over processors EMPI=1,EMPID


  !     write information on load-balancing to formatted file 'balance'

  !     INQUIRE(FILE='STOP',EXIST=STOPIT)
!  if (ITER.eq.SCFSTEPS) then
!
!    if (MYACTVRANK.eq.0) then
!
!      open (50,file='ebalance',form='formatted')
!      write(50, *) "# Energy load-balancing file"
!      write(50, *) "# 1st line: number of E-points, number of E-processes."
!      write(50, *) "# point --- process --- timing used"
!      write(50, *) IERLAST,EMPID
!
!      do IER = 1,IERLAST
!        write(50,*) IER, EPROC(IER), MTIME(IER)
!      end do
!
!      close(50)
!    endif
!
!  endif

end subroutine
