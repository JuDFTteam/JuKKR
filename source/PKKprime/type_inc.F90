!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


module type_inc

  implicit none

    type :: inc_TYPE

      sequence

      integer :: N = 33

      integer :: lmaxd  = -1
      integer :: lmax   = -1
      integer :: nspd   = -1
      integer :: nspo   = -1
      integer :: nspoh  = -1
      integer :: nrd    = -1
      integer :: nembd  = -1
      integer :: nemb   = -1
      integer :: nclsd  = -1
      integer :: ncls   = -1
      integer :: natypd = -1
      integer :: natyp  = -1
      integer :: naezd  = -1
      integer :: naez   = -1
      integer :: naclsd = -1
      integer :: ielast = -1
      integer :: ins    = -1

      integer :: lmmax  = -1
      integer :: lmmaxso= -1
      integer :: alm    = -1
      integer :: almso  = -1

      integer :: ndegen = -1
      integer :: nBZdim = -1
      integer :: nrootmax = -1

      integer          :: neig   = -1
      double precision :: reig   = 0.0
      logical          :: memopt = .false.
      logical          :: feast  = .false.

      logical :: lrhod = .false.
      logical :: ltorq = .false.
      logical :: lspinflux = .false.
      logical :: lalpha = .false.
!      logical :: simpson = .false.

    end type inc_TYPE

#ifdef CPP_MPI

contains

  subroutine create_mpimask_inc(inc,myMPItype,iout)

    use mpi
    implicit none

    type(inc_type), intent(in)  :: inc
    integer,        intent(out) :: myMPItype, iout

    integer :: N, blocklen(inc%N), etype(inc%N), ierr
    integer(kind=MPI_ADDRESS_KIND) :: disp(inc%N), base

    N = inc%N

    call MPI_Get_address(inc%N,      disp(1), ierr)
    call MPI_Get_address(inc%lmaxd,  disp(2), ierr)
    call MPI_Get_address(inc%lmax,   disp(3), ierr)
    call MPI_Get_address(inc%nspd,   disp(4), ierr)
    call MPI_Get_address(inc%nspo,   disp(5), ierr)
    call MPI_Get_address(inc%nspoh,  disp(6), ierr)
    call MPI_Get_address(inc%nrd,    disp(7), ierr)
    call MPI_Get_address(inc%nembd,  disp(8), ierr)
    call MPI_Get_address(inc%nemb,   disp(9), ierr)
    call MPI_Get_address(inc%nclsd,  disp(10), ierr)
    call MPI_Get_address(inc%ncls,   disp(11), ierr)
    call MPI_Get_address(inc%natypd, disp(12), ierr)
    call MPI_Get_address(inc%natyp,  disp(13), ierr)
    call MPI_Get_address(inc%naezd,  disp(14), ierr)
    call MPI_Get_address(inc%naez,   disp(15), ierr)
    call MPI_Get_address(inc%naclsd, disp(16), ierr)
    call MPI_Get_address(inc%ielast, disp(17), ierr)
    call MPI_Get_address(inc%ins,    disp(18), ierr)
    call MPI_Get_address(inc%lmmax,  disp(19), ierr)
    call MPI_Get_address(inc%lmmaxso,disp(20), ierr)
    call MPI_Get_address(inc%alm,    disp(21), ierr)
    call MPI_Get_address(inc%almso,  disp(22), ierr)
    call MPI_Get_address(inc%ndegen, disp(23), ierr)
    call MPI_Get_address(inc%nBZdim, disp(24), ierr)
    call MPI_Get_address(inc%nrootmax,disp(25), ierr)
    call MPI_Get_address(inc%neig,   disp(26), ierr)
    call MPI_Get_address(inc%reig,   disp(27), ierr)
    call MPI_Get_address(inc%memopt, disp(28), ierr)
    call MPI_Get_address(inc%feast,  disp(29), ierr)
    call MPI_Get_address(inc%lrhod,  disp(30), ierr)
    call MPI_Get_address(inc%ltorq,  disp(31), ierr)
    call MPI_Get_address(inc%lspinflux,  disp(32), ierr)
    call MPI_Get_address(inc%lalpha,  disp(33), ierr)
!    call MPI_Get_address(inc%simpson,  disp(28), ierr)

    base = disp(1)
    disp = disp - base

    blocklen=1

    etype(1:26)  = MPI_INTEGER
    etype(27)    = MPI_DOUBLE_PRECISION
    etype(28:33) = MPI_LOGICAL

    call MPI_Type_create_struct(N, blocklen, disp, etype, myMPItype, iout)

  end subroutine create_mpimask_inc
#endif


end module type_inc
