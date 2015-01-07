module mod_iohelp

  implicit none

  private
  public :: file_present, getBZvolume, get_nLines
#ifdef CPP_MPI
  public :: open_mpifile_setview, close_mpifile
#endif

contains

  logical function file_present(filename)

    use mod_mympi, only: myrank, nranks, master
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    character(len=*) :: filename
    logical          :: file_exist
    integer          :: ierr

    if(myrank==master) inquire(file=filename, exist=file_exist)

#ifdef CPP_MPI
    call MPI_Bcast( file_exist, 1, MPI_LOGICAL, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in Bcast(file_exist)'
#endif

    file_present = file_exist

  end function file_present



  function get_nLines( UIO ) result( nLines )
    implicit none
    integer :: nLines, IOERR, UIO
    character(len=256) :: tmpstring
    nLines = 0
    Read_Loop: DO
      READ( UIO, *, IOSTAT = IOERR ) tmpstring
      IF (IOERR < 0) THEN
        EXIT Read_Loop
      END IF
      nLines = nLines + 1
    END DO Read_Loop
  end function get_nLines


#ifdef CPP_MPI
  subroutine open_mpifile_setview(fileprefix, rwmode, ndim, dimens, nkpts, my_mpi_datatype, comm_myrank, comm_nranks, my_mpi_comm, filehandle)
    use mpi
    use mod_ioformat, only: fmt_fn_ext, ext_mpiio
    implicit none

    character(len=*), intent(in) :: fileprefix, rwmode
    integer, intent(in)  :: ndim, my_mpi_datatype, comm_myrank, comm_nranks, my_mpi_comm
    integer, intent(in)  :: dimens(ndim), nkpts(0:comm_nranks-1)
    integer, intent(out) :: filehandle

    integer :: myMPI_iotype
    integer(kind=MPI_OFFSET_KIND) :: disp=0, datasize=0
    integer, parameter :: mo = kind(disp)

    integer :: ihelp, nkpt, ntot, ioff, irank, ierr
    character(len=256) :: filename

    !input dependent parameters
    select case(my_mpi_datatype)
    case(MPI_INTEGER);          datasize=4_mo
    case(MPI_REAL);             datasize=4_mo
    case(MPI_COMPLEX);          datasize=8_mo
    case(MPI_DOUBLE_PRECISION); datasize=8_mo
    case(MPI_DOUBLE_COMPLEX);   datasize=16_mo
    case default; stop 'Datasize for my_mpi_datatype not known'
    end select

    !init
    ihelp =  product(dimens)
    nkpt = nkpts(comm_myrank)
    ntot = sum(nkpts)
    ioff = 0
    do irank=0,comm_myrank-1
      ioff = ioff + nkpts(irank)
    end do

    !open the file
    write(filename,fmt_fn_ext) trim(fileprefix), ext_mpiio
    select case (rwmode)
    case ('read')
      call MPI_File_open( my_mpi_comm, trim(filename), MPI_MODE_RDONLY, &
                        & MPI_INFO_NULL, filehandle, ierr                  )
      if(ierr/=MPI_SUCCESS) stop 'error: MPI_File_open readmode'
    case ('write')
      call MPI_File_open( my_mpi_comm, trim(filename), MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                        & MPI_INFO_NULL, filehandle, ierr                                  )
      if(ierr/=MPI_SUCCESS) stop 'error: MPI_File_open writemode'
    case default; stop 'dont know how to open the mpi-file. only mode = read/write'
    end select

    !prepare the file view (see only the part of this task)
    call MPI_Type_Contiguous( ihelp*nkpt, my_mpi_datatype, myMPI_iotype, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error: MPI_Type in open_mpifile_setview'

    call MPI_Type_commit(myMPI_iotype,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error: MPI_Type_commit'

    disp = ihelp*(ioff*datasize)!byte
    call MPI_File_set_view( filehandle, disp, my_mpi_datatype,             &
                          & myMPI_iotype, 'native', MPI_INFO_NULL, ierr )
    if(ierr/=MPI_SUCCESS) stop 'error: MPI_File_set_view'

  end subroutine open_mpifile_setview





  subroutine close_mpifile(filehandle)
    use mpi
    implicit none
    integer, intent(inout) :: filehandle
    integer :: ierr

    call MPI_File_close(filehandle,ierr)

  end subroutine close_mpifile

#endif




  double precision function getBZvolume(recbv)
    use mod_mathtools, only: crossprod
    implicit none
    double precision, intent(in) :: recbv(3,3)
    double precision :: cross(3)

    call crossprod(recbv(:,1), recbv(:,2), cross)
    getBZvolume = sum(cross*recbv(:,3))

  end function getBZvolume



end module mod_iohelp
