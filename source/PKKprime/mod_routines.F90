!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


module mod_routines

  implicit none

  private
  public :: routines

contains

  subroutine routines()

    use type_inc,  only: inc_type
    use type_data, only: lattice_type, cluster_type, tgmatrx_type
    use mod_mympi, only: mympi_init, myrank, nranks, master
    use mod_read,  only: read_inc, read_TBkkrdata
    use mod_fermisurf_basic,  only: testpath
    use mod_fermisurf,  only: fermisurface
    use mod_calconfs,   only: calc_on_fsurf_inputcard
    use mod_scattering, only: calc_scattering_inputcard

#ifdef CPP_MPI
    use mpi
#endif
#ifdef CPP_TIMING
    use mod_timing,     only: timing_init, timing_start, timing_stop
    use mod_types, only: t_inc ! needed to set i_time=0 for all but the master rank
#endif

    implicit none

    type(inc_type)      :: inc
    type(lattice_type)  :: lattice
    type(cluster_type)  :: cluster
    type(tgmatrx_type)  :: tgmatrx
    integer             :: ierr, nrootmax

    integer :: nkpts
    double precision, allocatable :: kpoints(:,:), areas(:)

    !initialize MPI
#ifdef CPP_MPI
    call MPI_Init ( ierr )
#endif
    call mympi_init()

#ifdef CPP_TIMING
    if (myrank/=master) t_inc%i_time = 0 ! diable writeing of the timing file for all but the master rank
    call timing_init(myrank, disable_serial_number=.true.)
#endif

    !Read in TBKKR-data
#ifdef CPP_TIMING
    call timing_start('Read TBkkr-data')
#endif
    call read_inc(inc)
    call read_TBkkrdata(inc, lattice, cluster, tgmatrx)

    !Perform tests
    if(myrank==master) call testpath(inc, lattice, cluster, tgmatrx)
#ifdef CPP_TIMING
    call timing_stop('Read TBkkr-data')
#endif

    !Calculate the Fermi Surface
    call fermisurface(inc, lattice, cluster, tgmatrx, nkpts, kpoints, areas)

    call calc_on_fsurf_inputcard(inc, lattice, cluster, tgmatrx, nkpts, kpoints)

    call calc_scattering_inputcard(inc, lattice, cluster, tgmatrx)

#ifdef CPP_MPI
    call MPI_Finalize(ierr)
#endif


  end subroutine routines

end module mod_routines
