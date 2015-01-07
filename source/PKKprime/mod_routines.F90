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

    !Read in TBKKR-data
    call read_inc(inc)

    call read_TBkkrdata(inc, lattice, cluster, tgmatrx)

    !Perform tests
    if(myrank==master) call testpath(inc, lattice, cluster, tgmatrx)

    !Calculate the Fermi Surface
    call fermisurface(inc, lattice, cluster, tgmatrx, nkpts, kpoints, areas)

    call calc_on_fsurf_inputcard(inc, lattice, cluster, tgmatrx, nkpts, kpoints)

    call calc_scattering_inputcard(inc, lattice, cluster, tgmatrx)

#ifdef CPP_MPI
    call MPI_Finalize(ierr)
#endif


  end subroutine routines

end module mod_routines
