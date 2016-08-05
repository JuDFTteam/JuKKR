!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


module mod_parutils

  implicit none

  private

  public :: distribute_linear_on_tasks, parallel_quicksort, parallel_quicksort_int

contains

  subroutine distribute_linear_on_tasks(nranks, myrank, master, ntot, ntot_pT, ioff_pT, output)
    !distributes 'ntot' points on 'nranks' tasks and returns the number of points per task, ntot_pT, and the offsets of the tasks, ioff_pT.

    implicit none

    integer, intent(in)  :: nranks, myrank, master, ntot
    logical, intent(in)  :: output
    integer, intent(out) :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1)

    integer :: irest, irank

    ntot_pT = int(ntot/nranks)
    ioff_pT = int(ntot/nranks)*(/ (irank, irank=0,nranks-1) /)
    irest   = ntot-int(ntot/nranks)*nranks

    if(irest>0) then

      do irank=0,irest-1
        ntot_pT(irank) = ntot_pT(irank) + 1
        ioff_pT(irank) = ioff_pT(irank) + irank
      end do!irank

      do irank=irest,nranks-1
        ioff_pT(irank) = ioff_pT(irank) + irest
      end do!irank

    end if!irest>0

    if(myrank==master .and. output) then
      write(*,*) '=== DISTRIBUTION OF K-POINTS ON TASKS: ==='
      do irank=0,nranks-1
        write(*,'("Task ",I0," treats points ",I0," to ",I0,", #of points= ",I0)') irank, ioff_pT(irank)+1, ioff_pT(irank)+ntot_pT(irank), ntot_pT(irank)
      end do!irank
      write(*,*) '=========================================='
    end if!myrank==master

  end subroutine distribute_linear_on_tasks


  subroutine parallel_quicksort(nranks, myrank, master, ntot, values, sortindex)
    use mod_mathtools, only: bubblesort
    use mpi
    implicit none

    integer :: nranks, myrank, master, ntot
    double precision, intent(in) :: values(ntot)
    integer, intent(out) :: sortindex(ntot)

    integer :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1), recvcount(0:nranks-1), displs(0:nranks-1)
    double precision :: minmax(2), minim, maxim, pivot(0:nranks)

    integer :: nsub, sublistpick(ntot)
    integer,          allocatable :: subsortind(:), sortindlocal(:)
    double precision, allocatable :: subvalues(:)

    integer :: ierr, irank, ii, ivp
    double precision :: mypivot_up, mypivot_dn
    double precision, parameter :: eps=1d-12

    !1) find minimal and maximal value
    minim=minval(values)
    maxim=maxval(values)
!   call distribute_linear_on_tasks(nranks, myrank, master, ntot, ntot_pT, ioff_pT, .false.)
!   npts = ntot_pT(myrank)
!   ioff = ioff_pT(myrank)
!   minmax(1) = minval(values(ioff+1:ioff+npts))
!   minmax(2) = maxval(values(ioff+1:ioff+npts))
!   !programm the MPI-Reduce here!!!

    !2) creare pivot points
    pivot(0)=minim-eps
    do ii=1,nranks
      pivot(ii) = (maxim-minim)/nranks*ii+minim
    end do!ii
    pivot(nranks) = pivot(nranks)+eps

    mypivot_dn = pivot(myrank)
    mypivot_up = pivot(myrank+1)

    !3) create (task-local) sublist
    nsub=0
    do ii=1,ntot
      if(values(ii)>mypivot_dn .and. values(ii)<=mypivot_up)then
        nsub = nsub+1
        sublistpick(nsub) = ii
      end if
    end do

    allocate(subvalues(nsub), subsortind(nsub), sortindlocal(nsub), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating subvalues etc.'
    subvalues(:) = values(sublistpick(1:nsub))


    !4) sort task-local sublist
    call bubblesort(nsub, subvalues, subsortind)
    sortindlocal = sublistpick(subsortind)

    !5) gather information from the other tasks
    call MPI_Allgather(nsub,1,MPI_INTEGER,recvcount,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'Error in MPI_Allgather(nsub)'
    displs = 0
    do irank=1,nranks-1
      displs(irank) = displs(irank-1) + recvcount(irank-1)
    end do
    call MPI_Allgatherv(sortindlocal,nsub,MPI_INTEGER,sortindex,recvcount,displs,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'Error in MPI_Allgatherv(sortindlocal)'


    !testoutput
!   if(myrank==master) write(*,'("Pivot: ",18ES18.9)') pivot
!   write(*,'("myrank= ",I2,", nsub= ",I8)') myrank, nsub
!   do ii=1,ntot
!     write(500+myrank,'(I8,2ES25.16)') ii, values(ii), values(sortindex(ii))
!   end do
  end subroutine parallel_quicksort




  subroutine parallel_quicksort_int(nranks, myrank, master, ntot, values, sortindex)
    use mod_mathtools, only: bubblesort_int
    use mpi
    implicit none

    integer :: nranks, myrank, master, ntot
    integer, intent(in) :: values(ntot)
    integer, intent(out) :: sortindex(ntot)

    integer :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1), recvcount(0:nranks-1), displs(0:nranks-1)
    double precision :: minmax(2), minim, maxim, pivot(0:nranks)

    integer :: nsub, sublistpick(ntot)
    integer, allocatable :: subsortind(:), sortindlocal(:)
    integer, allocatable :: subvalues(:)

    integer :: ierr, irank, ii, ivp
    double precision :: mypivot_up, mypivot_dn
    double precision, parameter :: eps=1d-12

    !1) find minimal and maximal value
    minim=dble(minval(values))
    maxim=dble(maxval(values))
!   call distribute_linear_on_tasks(nranks, myrank, master, ntot, ntot_pT, ioff_pT, .false.)
!   npts = ntot_pT(myrank)
!   ioff = ioff_pT(myrank)
!   minmax(1) = minval(values(ioff+1:ioff+npts))
!   minmax(2) = maxval(values(ioff+1:ioff+npts))
!   !programm the MPI-Reduce here!!!

    !2) create pivot points
    pivot(0)=minim-eps
    do ii=1,nranks
      pivot(ii) = (maxim-minim)/nranks*ii+minim
    end do!ii
    pivot(nranks) = pivot(nranks)+eps

    mypivot_dn = pivot(myrank)
    mypivot_up = pivot(myrank+1)

    !3) create (task-local) sublist
    nsub=0
    do ii=1,ntot
      if(dble(values(ii))>mypivot_dn .and. dble(values(ii))<=mypivot_up)then
        nsub = nsub+1
        sublistpick(nsub) = ii
      end if
    end do

    allocate(subvalues(nsub), subsortind(nsub), sortindlocal(nsub), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating subvalues etc.'
    subvalues(:) = values(sublistpick(1:nsub))


    !4) sort task-local sublist
    call bubblesort_int(nsub, subvalues, subsortind)
    sortindlocal = sublistpick(subsortind)

    !5) gather information from the other tasks
    call MPI_Allgather(nsub,1,MPI_INTEGER,recvcount,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'Error in MPI_Allgather(nsub)'
    displs = 0
    do irank=1,nranks-1
      displs(irank) = displs(irank-1) + recvcount(irank-1)
    end do
    call MPI_Allgatherv(sortindlocal,nsub,MPI_INTEGER,sortindex,recvcount,displs,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'Error in MPI_Allgatherv(sortindlocal)'


!   !testoutput
!   if(myrank==master) write(*,'("Pivot: ",18ES18.9)') pivot
!   write(*,'("myrank= ",I2,", nsub= ",I8)') myrank, nsub
!   do ii=1,ntot
!     write(500+myrank,'(I8,2I8)') ii, values(ii), values(sortindex(ii))
!   end do
  end subroutine parallel_quicksort_int


end module mod_parutils
