module mod_fermisurf

  implicit none

  private
  public :: fermisurface

contains

  subroutine fermisurface(inc, lattice, cluster, tgmatrx, nkpts_int, kpoints_int, areas_int)

    use type_inc,       only: inc_type
    use type_data,      only: lattice_type, cluster_type, tgmatrx_type
    use mod_symmetries, only: symmetries_type, set_symmetries
    use mod_mympi,      only: myrank, nranks, master
    use mod_iohelp,     only: file_present
    use mod_fermisurf_3D,    only: find_fermisurface_3D
    use mod_fermisurf_2D,    only: find_fermisurface_2D
    use mod_fermisurf_basic, only: save_kpointsfile_int

    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(cluster_type), intent(in) :: cluster
    type(tgmatrx_type), intent(inout) :: tgmatrx

    integer,                       intent(out) :: nkpts_int
    double precision, allocatable, intent(out) :: kpoints_int(:,:), areas_int(:)


    integer :: lfsurf, nCub3(3), nFSiter, nROOTiter, lrefine, nrefinenew
    integer, allocatable :: nCut_iter(:), roottype(:), nstepsconnect(:)
    double precision, allocatable :: rooteps(:)

    type(symmetries_type) :: symmetries



    call read_fscfg(lfsurf, nCub3, nFSiter, nROOTiter, nstepsconnect, nCut_iter, roottype, rooteps, lrefine, nrefinenew)
    if(lfsurf/=1) return

    !initialize symmetries
    call set_symmetries(inc, lattice, symmetries)

    select case (inc%nBZdim)
      case(2);  call find_fermisurface_2D( inc, lattice, cluster, tgmatrx, symmetries, nCub3, nFSiter, nROOTiter, nstepsconnect, &
                                         & nCut_iter, roottype, rooteps, lrefine, nrefinenew, nkpts_int, kpoints_int, areas_int  )
      case(3);  call find_fermisurface_3D( inc, lattice, cluster, tgmatrx, symmetries, nCub3, nFSiter, nROOTiter, nstepsconnect, &
                                         & nCut_iter, roottype, rooteps, lrefine, nrefinenew, nkpts_int, kpoints_int, areas_int  )
      case default; stop 'dimens must be 2 or 3'
    end select


!   save the integration k-points and areas to a file
    call save_kpointsfile_int(nkpts_int, kpoints_int, areas_int, symmetries%nsym_used, symmetries%isym_used)

  end subroutine fermisurface







  subroutine read_fscfg(lfsurf,nCub3,nFSiter,nROOTiter,nstepsconnect,nCut_iter,roottype,rooteps,lrefine,nrefinenew)
  !++++++++++++++++++++++++++++++++++++++++++++
  !+ read the parameter for the calculation   +
  !+ of the fermi surface from the input file +
  !++++++++++++++++++++++++++++++++++++++++++++
    use mod_ioinput, only: IoInput
    use mod_mympi,   only: myrank, nranks, master
    use mod_fermisurf_basic, only: ROOT_IMAG, ROOT_ANY, ROOT_REAL
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    integer, intent(out) :: lfsurf, nCub3(3), nFSiter, nROOTiter, lrefine, nrefinenew
    integer, allocatable, intent(out) :: nCut_iter(:), roottype(:), nstepsconnect(:)
    double precision, allocatable, intent(out) :: rooteps(:)

    integer              :: itmp(8), ierr, ii
    character(len=80)    :: uio
    character(len=8)     :: rootinp

    if(myrank==master) then
      call IoInput('LFSURF    ',uio,1,7,ierr)
      read(unit=uio,fmt=*) lfsurf

      call IoInput('NKPTCUBES ',uio,1,7,ierr)
      read(unit=uio,fmt=*) nCub3(:)

      call IoInput('NFSITER   ',uio,1,7,ierr)
      read(unit=uio,fmt=*) nFSiter

      allocate(nCut_iter(nFSiter), roottype(nFSiter+1), rooteps(nFSiter+1), nstepsconnect(nFSiter+1), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating nCut_iter etc.'

      call IoInput('NREFINE   ',uio,1,7,ierr)
      read(unit=uio,fmt=*) nCut_iter
!     write(1000,*) nCut_iter

      call IoInput('NROOTITER ',uio,1,7,ierr)
      read(unit=uio,fmt=*) nROOTiter

      do ii=1,nFSiter+1
        call IoInput('ROOTTYPE  ',uio,ii,7,ierr)
        read (unit=uio,fmt=*) rootinp
        select case (trim(rootinp))
          case('any');  roottype(ii)=ROOT_ANY
          case('real'); roottype(ii)=ROOT_REAL
          case('imag'); roottype(ii)=ROOT_IMAG
          case default; stop 'case for roottype not known: any/real/imag'
        end select
      end do!ii

      do ii=1,nFSiter+1
        call IoInput('ROOTEPS   ',uio,ii,7,ierr)
        read (unit=uio,fmt=*) rooteps(ii)
      end do!ii

      do ii=1,nFSiter+1
        call IoInput('NSTEPCCON ',uio,ii,7,ierr)
        read (unit=uio,fmt=*) nstepsconnect(ii)
      end do!ii

      if(lfsurf==1)then
        call IoInput('LREFINE   ',uio,ii,7,ierr)
        read (unit=uio,fmt=*) lrefine
        if(lrefine==1)then
          call IoInput('NREFINE2  ',uio,ii,7,ierr)
          read (unit=uio,fmt=*) nrefinenew
!       write(*,'("NREFINE= ",I0)') nrefinenew
        end if!lrefine==1
      end if!lfsurf

    end if!myrank==master

#ifdef CPP_MPI
    itmp(1)   = lfsurf
    itmp(2:4) = nCub3
    itmp(5)   = nFSiter
    itmp(6)   = nROOTiter
    itmp(7)   = lrefine
    itmp(8)   = nrefinenew

    call MPI_Bcast(itmp,8,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem broadcasting itmp'

    lfsurf        = itmp(1)
    nCub3         = itmp(2:4)
    nFSiter       = itmp(5)
    nROOTiter     = itmp(6)
    lrefine       = itmp(7)
    nrefinenew    = itmp(8)

    if(myrank/=master) then
      allocate(nCut_iter(nFSiter), roottype(nFSiter+1), rooteps(nFSiter+1), nstepsconnect(nFSiter+1), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating nCut_iter etc.'
    end if!myrank/=master

    call MPI_Bcast(nCut_iter,nFSiter,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem broadcasting nCut_iter'

    call MPI_Bcast(roottype,nFSiter+1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem broadcasting nCut_iter'

    call MPI_Bcast(rooteps,nFSiter+1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem broadcasting nCut_iter'

    call MPI_Bcast(nstepsconnect,nFSiter+1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem broadcasting nCut_iter'

#endif

  end subroutine read_fscfg






end module mod_fermisurf
