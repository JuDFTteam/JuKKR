!--------------------------------------------------------------------------------
! Copyright (c) 2018 Forschungszentrum Juelich GmbH, Juelich, Germany
! This file is part of KKRnano and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module Statistics_mod
!-------------------------------------------------------------------------------
!> Summary: Reconstruction of a continous distribution function without prior knowledge of the range of data
!> Author: Paul F Baumeister
!> Category: KKRnano
!-------------------------------------------------------------------------------
  implicit none
  private
  
  public :: SimpleStats, init, add, allreduce, eval
  
  ! planned: getMPIstats --> each process adds his number. 
  !          Can already be emulated by this sequence: init(s); add(s, mynumber); i=allreduce(s, comm); eval(s)

  integer, public, parameter :: NUMBER_OF_MOMENTS = 4
  
  type SimpleStats
    real(kind=8) :: moment(0:NUMBER_OF_MOMENTS-1) ! [size(x), sum(x), sum(x**2), ...] ! MPI_Allreduce with MPI_SUM
    real(kind=8) :: xtrema(0:1) ! [-minval(x), maxval(x)], minus minval is stored, so MPI_Allreduce both with MPI_MAX 
    character(len=16) :: name
  endtype

  interface init
    module procedure init_stats
  endinterface
  
  interface add
    module procedure add_1r8, add_nr8
  endinterface

  interface allreduce
    module procedure allreduce_nss
  endinterface
  
  interface eval
    module procedure eval_stats
  endinterface
  
  contains

  elemental subroutine init_stats(self, name)
    type(SimpleStats), intent(inout) :: self
    character(len=*),  intent(in)    :: name
    self%moment = 0.d0 ! init
    self%xtrema = -huge(1.d0) ! init
    self%name = adjustl(name)
  endsubroutine ! init

  subroutine add_1r8(self, x)
    type(SimpleStats), intent(inout) :: self
    real(kind=8), intent(in) :: x

    self%moment(0:3) = self%moment(0:3) + [1.d0, x, x*x, x*x*x]
    self%xtrema(0:1) = max(self%xtrema(0:1), [-x, x])
  endsubroutine ! reduce

  subroutine add_nr8(self, x)
    type(SimpleStats), intent(inout) :: self(:)
    real(kind=8), intent(in) :: x(:)
    integer :: i ! group index
    do i = 1, size(self)
      self(i)%moment(0:3) = self(i)%moment(0:3) + [1.d0, x(i), x(i)*x(i), x(i)*x(i)*x(i)]
      self(i)%xtrema(0:1) = max(self(i)%xtrema(0:1), [-x(i), x(i)])
    enddo ! i
  endsubroutine ! reduce

  integer function allreduce_nss(self, communicator) result(ierr)
    include 'mpif.h'
    type(SimpleStats), intent(inout) :: self(:) ! group many to save communication latencies
    integer, intent(in) :: communicator

    real(kind=8), allocatable :: ssum(:,:,:), smax(:,:,:)
    integer, parameter :: iSEND=1, iRECV=2, m = NUMBER_OF_MOMENTS
    integer :: ist, ngroup, i

    ngroup = size(self)
    allocate(ssum(0:m-1,ngroup,iSEND:iRECV), smax(0:1,ngroup,iSEND:iRECV), stat=ierr)
    
    ! pack
    do i = 1, ngroup
      ssum(:,i,iSend) = self(i)%moment(:) 
      smax(:,i,iSend) = self(i)%xtrema(:) 
    enddo ! i
    
    call MPI_Allreduce(ssum(:,:,iSend), ssum(:,:,iRECV), m*ngroup, MPI_REAL8, MPI_SUM, communicator, ierr)
    call MPI_Allreduce(smax(:,:,iSend), smax(:,:,iRECV), 2*ngroup, MPI_REAL8, MPI_MAX, communicator, ierr)
    
    ! unpack
    do i = 1, ngroup
      self(i)%moment(:) = ssum(:,i,iRECV) 
      self(i)%xtrema(:) = smax(:,i,iRECV) 
    enddo ! i

    deallocate(ssum, smax, stat=ist) ! ignore status
  endfunction ! allreduce

  character(len=96) function eval_stats(self, scaleby, ndigits) result(str)
    type(SimpleStats),      intent(in) :: self
    real(kind=8), optional, intent(in) :: scaleby ! scale results by this
    integer,      optional, intent(in) :: ndigits ! if we know that the data are all integer, we should pass ndigits=0
    
    integer :: ios, ndig!, nd10
    character(len=32) :: frmt, minmax
    real(kind=8) :: invN, mean, var, dev, xtr(0:1), scal
    
    scal = 1.d0 ; if (present(scaleby)) scal = scaleby
    ndig = 1    ; if (present(ndigits)) ndig = ndigits

    invN = 1.d0/dble(max(1.d0, self%moment(0)))
    mean = self%moment(1)*invN ! divide by the number of samples
    
!   nd10 = ceiling(log10(mean))

    var  = 0 ; if (ubound(self%moment, 1) > 1) &
    var  = self%moment(2)*invN - mean*mean ! variance
    dev  = sqrt(max(0.d0, var)) ! compute sigma as sqrt(variance)
    xtr(0:1) = self%xtrema(0:1)*[-scal, scal]
    if (all(abs(xtr - nint(xtr)) < 0.1**(ndig+1))) then
      write(unit=minmax, fmt='("  [",i0,", ",i0,"]")', iostat=ios) nint(xtr)   
    else
      write(unit=frmt, fmt="(9(a,i0))", iostat=ios) '("  [",f0.",ndig,",", ",f0.",ndig,"]")' ! generate format string
      write(unit=minmax, fmt=frmt, iostat=ios) xtr
    endif
    write(unit=frmt, fmt="(9(a,i0))", iostat=ios) "(a,9(a,f0.",ndig,"))" ! generate format string
    write(unit=str, fmt=frmt, iostat=ios) trim(self%name),": ",mean*scal," +/- ",dev*scal, trim(minmax)
  endfunction ! eval
 
endmodule ! Statistics_mod
