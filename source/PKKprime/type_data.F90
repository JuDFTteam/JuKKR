!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


module type_data

  implicit none

    private
    public :: lattice_TYPE, lattice_init,&
            & cluster_TYPE, cluster_init,&
            & tgmatrx_TYPE, tgmatrx_init,&
            & get_lattice

#ifdef CPP_MPI
    public :: create_mpimask_lattice, create_mpimask_cluster, create_mpimask_tgmatrx
#endif

    type :: lattice_TYPE

      integer          :: N = 5

      double precision :: alat
      double precision :: bravais(3,3)
      double precision :: recbv(3,3)
      double precision, allocatable :: rbasis(:,:)

    end type lattice_TYPE


    type :: cluster_TYPE

      integer          :: N = 7

      integer,          allocatable :: cls(:)
      integer,          allocatable :: nacls(:)
      integer,          allocatable :: ezoa(:,:)
      integer,          allocatable :: atom(:,:)
      double precision, allocatable :: rcls(:,:,:)
      double precision, allocatable :: rr(:,:)

    end type cluster_TYPE


    type :: tgmatrx_TYPE

      integer          :: N = 10

      double complex, allocatable :: energies(:)
      double complex, allocatable :: tmatll(:,:,:,:)
      double complex, allocatable :: tmat(:,:,:)
      double complex, allocatable :: tinvll(:,:,:,:)
      double complex, allocatable :: ginp(:,:,:,:)
      double complex, allocatable :: rhod(:,:,:,:)
      double complex, allocatable :: torq(:,:,:,:)
      double complex, allocatable :: spinflux(:,:,:,:)
      double complex, allocatable :: alpha(:,:,:,:)

    end type tgmatrx_TYPE



contains



  subroutine lattice_init(inc,lattice)

    use type_inc, only: inc_TYPE
    implicit none

    type(inc_TYPE),     intent(in)  :: inc
    type(lattice_TYPE), intent(out) :: lattice

    integer :: ierr

    lattice%alat         = -1d0
    lattice%bravais(:,:) = -1d0
    lattice%recbv(:,:)   = -1d0

    allocate(lattice%rbasis(3,inc%naezd+inc%nembd), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating lattice%rbasis'

    lattice%rbasis(:,:) = -1d0

  end subroutine lattice_init



  function get_lattice( bravais ) result ( lattype )

    implicit none

    double precision, intent(in) :: bravais(3,3)
    character(len=3) :: lattype

    double precision :: fccdist, bccdist, hexdist
    double precision, parameter :: eps=1d-7
    double precision, parameter :: bcclat(3,3)=reshape( (/ -.5d0, .5d0, .5d0, .5d0, -.5d0, .5d0, .5d0, .5d0, -.5d0 /), (/ 3, 3/) ), &
                                 & fcclat(3,3)=reshape( (/   0d0, .5d0, .5d0, .5d0,   0d0, .5d0, .5d0, .5d0,   0d0 /), (/ 3, 3/) ), &
                                 & hexlat(2,2)=reshape( (/  .5d0, -0.8660254038d0, .5d0, 0.8660254038d0 /), (/ 2,2 /) )

    fccdist = sqrt(sum( (bravais-fcclat)**2 ))/9
    bccdist = sqrt(sum( (bravais-bcclat)**2 ))/9
    hexdist = sqrt(sum( (bravais(1:2,1:2)-hexlat)**2 ))/4

    if(fccdist<eps)then
      lattype = 'fcc'
    elseif(bccdist<eps)then
      lattype = 'bcc'
    elseif(hexdist<eps)then
      lattype = 'hex'
    else
      stop 'in get_lattice: lattice not recognized'
    end if

  end function get_lattice



  subroutine cluster_init(inc,cluster)

    use type_inc, only: inc_TYPE
    implicit none

    type(inc_TYPE),     intent(in)  :: inc
    type(cluster_TYPE), intent(out) :: cluster

    integer :: ierr

    allocate( cluster%cls(inc%natypd),              &
            & cluster%nacls(inc%nclsd),             &
            & cluster%ezoa(inc%naclsd,inc%naezd),   &
            & cluster%atom(inc%naclsd,inc%naezd),   &
            & cluster%rcls(3,inc%naclsd,inc%nclsd), &
            & cluster%rr(3,0:inc%nrd),              &
            & STAT=ierr                     )
    if(ierr/=0) stop 'Problem allocating arrays for cluster'

    cluster%cls   = -1
    cluster%nacls = -1
    cluster%ezoa  = -1
    cluster%atom  = -1
    cluster%rcls  = -1d0
    cluster%rr    = -1d0

  end subroutine cluster_init



  subroutine tgmatrx_init(inc,tgmatrx)

    use type_inc, only: inc_TYPE
    implicit none

    type(inc_TYPE),     intent(in)  :: inc
    type(tgmatrx_TYPE), intent(out) :: tgmatrx

    integer :: ierr

    allocate( tgmatrx%energies(inc%ielast),                                     &
            & tgmatrx%tmatll(inc%lmmaxso,inc%lmmaxso,inc%naezd,inc%ielast),     &
            & tgmatrx%tmat(inc%almso,inc%almso,inc%ielast),                     &
            & tgmatrx%tinvll(inc%lmmaxso,inc%lmmaxso,inc%naezd,inc%ielast),     &
            & tgmatrx%ginp(inc%naclsd*inc%lmmax,inc%lmmax,inc%nclsd,inc%ielast),&
            & tgmatrx%rhod(inc%lmmaxso,inc%lmmaxso,inc%natypd,4),               &
            & tgmatrx%torq(inc%lmmaxso,inc%lmmaxso,inc%natypd,3),               &
            & tgmatrx%spinflux(inc%lmmaxso,inc%lmmaxso,inc%natypd,3),           &
            & tgmatrx%alpha(inc%lmmaxso,inc%lmmaxso,inc%natypd,3),              &
            & STAT=ierr                                                         )
    if(ierr/=0) stop 'Problem allocating arrays for tgmatrx'

    tgmatrx%energies = (-1d0,0d0)
    tgmatrx%tmatll   = (-1d0,0d0)
    tgmatrx%tmat     = (-1d0,0d0)
    tgmatrx%tinvll   = (-1d0,0d0)
    tgmatrx%ginp     = (-1d0,0d0)
    tgmatrx%rhod     = (0d0,0d0)
    tgmatrx%torq     = (0d0,0d0)
    tgmatrx%spinflux = (0d0,0d0)
    tgmatrx%alpha    = (0d0,0d0)

  end subroutine tgmatrx_init


#ifdef CPP_MPI

  subroutine create_mpimask_lattice(lattice,myMPItype)

    use mpi
    implicit none

    type(lattice_type), intent(in)  :: lattice
    integer,            intent(out) :: myMPItype

    integer :: blocklen(lattice%N), etype(lattice%N), ierr
    integer(kind=MPI_ADDRESS_KIND) :: disp(lattice%N), base

    call MPI_Get_address(lattice%N,       disp(1), ierr)
    call MPI_Get_address(lattice%alat,    disp(2), ierr)
    call MPI_Get_address(lattice%bravais, disp(3), ierr)
    call MPI_Get_address(lattice%recbv,   disp(4), ierr)
    call MPI_Get_address(lattice%rbasis,  disp(5), ierr)

    base = disp(1)
    disp = disp - base

    blocklen(1)=1
    blocklen(2)=1
    blocklen(3)=size(lattice%bravais)
    blocklen(4)=size(lattice%recbv)
    blocklen(5)=size(lattice%rbasis)

    etype(1)   = MPI_INTEGER
    etype(2:5) = MPI_DOUBLE_PRECISION

    call MPI_Type_create_struct(lattice%N, blocklen, disp, etype, myMPItype, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in create_mpimask_lattice'

  end subroutine create_mpimask_lattice

  subroutine create_mpimask_cluster(cluster,myMPItype)

    use mpi
    implicit none

    type(cluster_type), intent(in)  :: cluster
    integer,            intent(out) :: myMPItype

    integer :: blocklen(cluster%N), etype(cluster%N), ierr
    integer(kind=MPI_ADDRESS_KIND) :: disp(cluster%N), base

    call MPI_Get_address(cluster%N,       disp(1), ierr)
    call MPI_Get_address(cluster%cls,     disp(2), ierr)
    call MPI_Get_address(cluster%nacls,   disp(3), ierr)
    call MPI_Get_address(cluster%ezoa,    disp(4), ierr)
    call MPI_Get_address(cluster%atom,    disp(5), ierr)
    call MPI_Get_address(cluster%rcls,    disp(6), ierr)
    call MPI_Get_address(cluster%rr,      disp(7), ierr)

    base = disp(1)
    disp = disp - base

    blocklen(1)=1
    blocklen(2)=size(cluster%cls)
    blocklen(3)=size(cluster%nacls)
    blocklen(4)=size(cluster%ezoa)
    blocklen(5)=size(cluster%atom)
    blocklen(6)=size(cluster%rcls)
    blocklen(7)=size(cluster%rr)

    etype(1:5)   = MPI_INTEGER
    etype(6:7) = MPI_DOUBLE_PRECISION

    call MPI_Type_create_struct(cluster%N, blocklen, disp, etype, myMPItype, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in create_mpimask_cluster'

  end subroutine create_mpimask_cluster


  subroutine create_mpimask_tgmatrx(tgmatrx,myMPItype)

    use mpi
    implicit none

    type(tgmatrx_type), intent(in)  :: tgmatrx
    integer,            intent(out) :: myMPItype

    integer :: blocklen(tgmatrx%N), etype(tgmatrx%N), ierr
    integer(kind=MPI_ADDRESS_KIND) :: disp(tgmatrx%N), base

    call MPI_Get_address(tgmatrx%N,       disp(1), ierr)
    call MPI_Get_address(tgmatrx%energies,disp(2), ierr)
    call MPI_Get_address(tgmatrx%tmatll,  disp(3), ierr)
    call MPI_Get_address(tgmatrx%tmat,    disp(4), ierr)
    call MPI_Get_address(tgmatrx%tinvll,  disp(5), ierr)
    call MPI_Get_address(tgmatrx%ginp,    disp(6), ierr)
    call MPI_Get_address(tgmatrx%rhod,    disp(7), ierr)
    call MPI_Get_address(tgmatrx%torq,    disp(8), ierr)
    call MPI_Get_address(tgmatrx%spinflux,disp(9), ierr)
    call MPI_Get_address(tgmatrx%alpha,   disp(10), ierr)

    base = disp(1)
    disp = disp - base

    blocklen(1)=1
    blocklen(2)=size(tgmatrx%energies)
    blocklen(3)=size(tgmatrx%tmatll)
    blocklen(4)=size(tgmatrx%tmat)
    blocklen(5)=size(tgmatrx%tinvll)
    blocklen(6)=size(tgmatrx%ginp)
    blocklen(7)=size(tgmatrx%rhod)
    blocklen(8)=size(tgmatrx%torq)
    blocklen(9)=size(tgmatrx%spinflux)
    blocklen(10)=size(tgmatrx%alpha)

    etype(1)   = MPI_INTEGER
    etype(2:10) = MPI_DOUBLE_COMPLEX

    call MPI_Type_create_struct(tgmatrx%N, blocklen, disp, etype, myMPItype, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in create_mpimask_tgmatrx'

  end subroutine create_mpimask_tgmatrx

#endif


end module type_data
