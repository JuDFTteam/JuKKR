!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


!-------------------------------------------------------------------------------
!> Summary: Module containing parallelization infos
!> Author: B. Zimmermann
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!-------------------------------------------------------------------------------
module mod_mympi

implicit none

  private
  public :: myrank, nranks, master, mympi_init

  integer, save :: myrank = -1
  integer, save :: nranks = -1
  integer, save :: master = -1

contains

  !-------------------------------------------------------------------------------
  !> Summary: Initialize MPI and set myrank, nranks
  !> Author: B. Zimmermann
  !> Category: PKKprime, communication
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Initializes myrank and nranks even if MPI is not used
  !-------------------------------------------------------------------------------
  subroutine mympi_init()

#ifdef CPP_MPI
    use mpi
#endif

    integer :: ierr

    master = 0

#ifdef CPP_MPI
    call MPI_Comm_rank ( MPI_COMM_WORLD, myrank, ierr )
    call MPI_Comm_size ( MPI_COMM_WORLD, nranks, ierr )
#else
    myrank = master
    nranks = 1
#endif

  end subroutine mympi_init


end module mod_mympi
