! TODO: problem:
!       the Fourier-transform needs all the ref. cluster data see kloopz1
! For building the coefficient matrix NUMN0 and INDN0 (off all atoms in LIZ)
! are needed!!!!!

! TODO: How to deal with lattice vectors???

! TODO: Separate into unique and non-unique atom clusters

module RefClusterData_mod
  implicit none

  type RefClusterData
    integer :: cluster_index
    !> index of lattice vector of non-unique cluster atoms
    integer, dimension(:), allocatable :: EZOA
    !> number of unique cluster atoms
    integer :: NUMN0
    !> atom indices of UNIQUE cluster atoms
    integer, dimension(:), allocatable :: INDN0
    !> atom indices of non-unique cluster atoms
    integer, dimension(:), allocatable :: ATOM
    !> number of non-unique cluster atoms
    integer :: naclsd
    !> positions of the non-unique cluster atoms relative to central atom
    double precision, dimension(:,:), allocatable :: RCLS
  end type

CONTAINS

  !----------------------------------------------------------------------------
  subroutine createRefClusterData(cluster, naclsd)
    implicit none
    type (RefClusterData), intent(inout) :: cluster
    integer, intent(in) :: naclsd

    cluster%naclsd = naclsd


    ! TODO: check numn0 <= naclsd

    ! TODO: OVERDIMENSIONED! SHOULD be numn0!!!
    ! But: use naclsd to keep combatibility to old prog.
    allocate(cluster%indn0(naclsd))

    allocate(cluster%ezoa(naclsd))
    allocate(cluster%atom(naclsd))
    allocate(cluster%RCLS(3, naclsd))

    cluster%cluster_index = -1
    cluster%indn0 = -1
    cluster%numn0 = -1
    cluster%ezoa = -1
    cluster%atom = -1
    cluster%RCLS = 0.0d0

  end subroutine

  !----------------------------------------------------------------------------
  subroutine destroyRefClusterData(cluster)
    implicit none
    type (RefClusterData), intent(inout) :: cluster

    deallocate(cluster%indn0)
    deallocate(cluster%ezoa)
    deallocate(cluster%atom)
    deallocate(cluster%RCLS)

  end subroutine

  !----------------------------------------------------------------------------
  !> Write reference cluster data to direct access file 'fileunit'
  !> at record 'recnr'
  subroutine writeRefClusterDataDA(cluster, fileunit, recnr)

    implicit none
    type (RefClusterData), intent(in) :: cluster
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr

    integer, parameter :: MAGIC_NUMBER = -1342263298

    write (fileunit, rec=recnr) MAGIC_NUMBER, &
                                cluster%cluster_index, &
                                cluster%indn0, &
                                cluster%numn0, &
                                cluster%ezoa, &
                                cluster%atom, &
                                cluster%RCLS, &
                                MAGIC_NUMBER
  end subroutine

  !----------------------------------------------------------------------------
  !> Read ref. cluster data from direct access file 'fileunit' at record 'recnr'
  subroutine readRefClusterDataDA(cluster, fileunit, recnr)
    implicit none

    type (RefClusterData), intent(inout) :: cluster
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr

    integer, parameter :: MAGIC_NUMBER = -1342263298
    integer :: magic, magic2

    read  (fileunit, rec=recnr) magic, &
                                cluster%cluster_index, &
                                cluster%indn0, &
                                cluster%numn0, &
                                cluster%ezoa, &
                                cluster%atom, &
                                cluster%RCLS, &
                                magic2

    if (magic /= MAGIC_NUMBER .or. magic2 /= MAGIC_NUMBER) then
      write (*,*) "ERROR: Invalid cluster data read. ", __FILE__, __LINE__
      STOP
    end if
  end subroutine

  !----------------------------------------------------------------------------
  !> Opens RefClusterData direct access file.
  subroutine openRefClusterDataDAFile(cluster, fileunit, filename)
    implicit none

    type (RefClusterData), intent(in) :: cluster
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename

    !------
    integer :: reclen

    integer, parameter :: MAGIC_NUMBER = -1342263298

    inquire (iolength = reclen) MAGIC_NUMBER, &
                                cluster%cluster_index, &
                                cluster%indn0, &
                                cluster%numn0, &
                                cluster%ezoa, &
                                cluster%atom, &
                                cluster%RCLS, &
                                MAGIC_NUMBER

    open(fileunit, access='direct', file=filename, recl=reclen, form='unformatted')

  end subroutine

  !----------------------------------------------------------------------------
  !> Closes RefClusterData direct access file.
  subroutine closeRefClusterDataDAFile(fileunit)
    implicit none
    integer, intent(in) :: fileunit

    close(fileunit)

  end subroutine


end module
