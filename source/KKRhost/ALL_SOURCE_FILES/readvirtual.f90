module mod_readvirtual

contains

program readvirtual
  use :: mod_datatypes, only: dp
  implicit none
  integer, parameter :: naezd = 20
  integer :: naez, naeznew
  integer :: natomv
  real (kind=dp), allocatable :: ratomv(:, :)
  real (kind=dp) :: rbasislist(3, naezd)
  real (kind=dp) :: rbasisold(3, naezd)
  integer :: iatom, ibasis
  integer :: ierr
  real (kind=dp) :: ratomvtest(3)
  real (kind=dp) :: rbasis(3), rbasisnew(3)

  real (kind=dp) :: bravais(3, 3)
  integer :: ndim

  ! get Bravais vectors <- implement !

  ! read the scoef file
  open (unit=32452345, file='scoef', iostat=ierr)
  if (ierr/=0) stop '[readvirtual] file not found'
  read (32452345, *) natomv
  allocate (ratomv(3,natomv))
  do iatom = 1, natomv
    read (32452345, *) ratomv(:, iatom)
  end do


  ! IF ( LSURF ) THEN
  ! NDIM = 2
  ! WRITE (6,'(23X,A)') 'LATTIX99: surface geometry mode'
  ! ELSE
  ! NDIM = 3
  ! WRITE (6,'(23X,A)') '  LATTIX99: bulk geometry mode'
  ! END IF
  ! DO I = 1,NDIM
  ! CALL IOINPUT('BRAVAIS   ',UIO,I,7,IER)
  ! READ (UNIT=UIO,FMT=*) (BRAVAIS(J,I),J=1,NDIM)
  ! END DO



  bravais(1, :) = (/ 1.0d0, 0.0d0, 0.0d0 /)
  bravais(2, :) = (/ 0.0d0, 1.0d0, 0.0d0 /)
  bravais(3, :) = (/ 0.0d0, 0.0d0, 1.0d0 /)


  rbasisold(:, 1) = (/ -0.9d0, -0.9d0, -0.9d0 /)
  rbasisold(:, 2) = (/ 0.15d0, 0.15d0, 0.15d0 /)
  rbasisold(:, 3) = (/ 0.2d0, 0.2d0, 0.2d0 /)
  naez = 3



  if (bravais(1,3)**2+bravais(2,3)**2+bravais(3,3)**2<10e-10) then
    ndim = 2
  else
    ndim = 3
  end if

  do ibasis = 1, naez
    call rtobasis(bravais, rbasisold(:,ibasis), rbasisnew, ndim)
    rbasislist(:, ibasis) = rbasisnew
  end do                           ! iatom
  ibasis = naez + 1


  do iatom = 1, natomv
    call rtobasis(bravais, ratomv(:,iatom), rbasisnew, ndim)
    if (.not. vec_in_list(rbasisnew,rbasislist,ibasis)) then
      rbasislist(:, ibasis) = rbasisnew
      ibasis = ibasis + 1
    end if
  end do                           ! natomv
  naeznew = ibasis - 1

  do ibasis = 1, naeznew
    write (*, *) rbasislist(:, ibasis)
  end do
  ! ratomvtest = (/ 1.5D0,1.5D0,1.5D0/)

  ! call rtobasis(bravais,ratomvtest,rbasis,ndim)
  ! write(*,*) 'RBASIS is ',rbasis


contains

  logical function vec_in_list(vec, veclist, bound)
    ! --------------------------
    ! checks if the vector vec is in the vector list veclist
    ! in the range of (1,bound)
    ! --------------------------
    integer :: bound
    real (kind=dp) :: vec(3)
    real (kind=dp) :: veclist(3, bound)
    integer :: ilist
    real (kind=dp) :: tempvec(3), diff

    vec_in_list = .false.
    do ilist = 1, bound
      tempvec = vec - veclist(:, ilist)
      diff = sqrt(tempvec(1)**2+tempvec(2)**2+tempvec(3)**2)
      if (diff<10e-10) vec_in_list = .true.
    end do
  end function vec_in_list


  subroutine rtobasis(bravais, rpos, rbasis, ndim)
    ! --------------------------
    ! converts a spacial vector rpos to a basis vector rbasis
    ! such that rbasis = bravais * n with n in [0,1]^ndim
    ! --------------------------
    implicit none
    real (kind=dp) :: bravais(3, 3)
    real (kind=dp) :: bravais_inv(ndim, ndim)
    integer :: ndim
    real (kind=dp) :: rpos(3)
    real (kind=dp) :: rbasis(3)

    real (kind=dp) :: ncoeffreal(ndim)
    integer :: ncoeffint(ndim)
    integer :: idim

    ! --------------------------
    ! first invert the bravais matrix => bravais_inv
    ! --------------------------
    bravais_inv = bravais(1:ndim, 1:ndim)
    call inverse_d1(bravais_inv)
    ! --------------------------
    ! then do n = Bravais* rpos
    ! --------------------------
    ncoeffreal = 0.0d0
    do idim = 1, ndim
      ncoeffreal = ncoeffreal + bravais_inv(1:ndim, idim)*rpos(idim)
    end do
    ! --------------------------
    ! take the smaller integer value of n => ncoeffint
    ! --------------------------
    do idim = 1, ndim
      ncoeffint(idim) = floor(ncoeffreal(idim))
    end do
    ! --------------------------
    ! do rbasis = rpos - bravias * ncoeffint
    ! --------------------------
    rbasis = 0.0d0
    do idim = 1, ndim
      rbasis = rbasis + bravais(:, idim)*ncoeffint(idim)
    end do
    rbasis = rpos - rbasis
  end subroutine rtobasis

  subroutine inverse_d1(mat)
    implicit none
    real (kind=8), intent (inout) :: mat(:, :)
    real (kind=8), allocatable :: work(:)
    integer, allocatable :: ipiv(:)
    integer :: n, info, i, j

    n = size(mat, 1)
    if (size(mat,2)/=n) stop 'inverse_d1: array dimensions differ.'
    allocate (ipiv(n), work(n))
    call dgetrf(n, n, mat, n, ipiv, info)
    if (info/=0) stop 'inverse_d1: dpotrf failed.'
    call dgetri(n, mat, n, ipiv, work, n, info)
    if (info/=0) stop 'inverse_d1: dpotri failed.'
  end subroutine inverse_d1


end program readvirtual

end module mod_readvirtual
