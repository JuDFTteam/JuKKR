!------------------------------------------------------------------------------------
!> Summary: Provides functionality to add virtual sites to the lattice
!> Author: 
!------------------------------------------------------------------------------------
module mod_addvirtual14

  ! by default everything is private
  private

  public :: addviratoms14

contains

  !-------------------------------------------------------------------------------
  !> Summary: Adds virtual atom sites to basis
  !> Author: 
  !> Category: KKRhost, geometry, input-output
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !> 
  !> Reads position of virtual atoms (used by impurity code to place atoms there) from scoef file 
  !> and adds virtual positions to basis. Automatically sets corresponding refpot etc. values
  !-------------------------------------------------------------------------------
  subroutine addviratoms14(linterface, nvirt, naez, naezd, natypd, nemb, nembd, rbasis, lcartesian, bravais, ncls, nineq, refpot, kaoez, noq, nref, rmtrefat, i25)

    use :: mod_datatypes, only: dp
    use :: mod_getclusnxyz, only: getclusnxyz
    implicit none
    ! interface variables
    logical :: linterface, lcartesian, labscord
    integer :: naez, ncls, nineq
    integer :: naezd, nembd
    integer :: natypd, nref, nvirt
    integer :: nemb
    integer :: ivir(1000)
    integer, parameter :: naclsd=1000

    integer :: i, i1, j
    real (kind=dp) :: rbasis(3, *), rbasisold(3, nemb+naezd), rbasissave(3, nemb+naezd)
    real (kind=dp) :: rmtrefat(naezd+nembd)
    integer :: refpot(*), refpotold(naezd+nemb)
    integer :: noq(*)
    integer :: kaoez(natypd, *), kaoezold(1, nemb+naezd)
    real (kind=dp) :: diff, rmaxclus, vec1(3), vec2(3, naclsd)
    integer :: nbr(3), nmax, nmaxz, n1, n2, n3, iq

    ! local variables
    character (len=40) :: i25
    integer :: nrefold
    integer :: natomimp
    real (kind=dp), allocatable :: ratomimp(:, :)
    ! real (kind=dp),allocatable  :: rbasislist(:,:)
    integer, allocatable :: atomimp(:)
    integer :: iatom, ibasis
    integer :: ierr
    ! real (kind=dp)              :: ratomvtest(3)
    real (kind=dp) :: rbasisnew(3)
    real (kind=dp), allocatable :: rclsnew(:, :), rbasisnew1(:, :)

    real (kind=dp) :: bravais(3, 3)
    real (kind=dp), allocatable :: bravaisinv(:, :)
    real (kind=dp) :: tol
    integer :: ndim
    integer :: naeznew

    tol = 1.e-5_dp

    write (1337, *) 'LINTERFACE', linterface
    write (1337, *) 'NAEZ', naez
    write (1337, *) 'NAEZD', naezd
    write (1337, *) 'NEMB', nemb
    write (1337, *) 'RBASISOLD'
    do ibasis = 1, naez + nemb
      write (1337, *) ibasis, rbasis(:, ibasis)
      refpotold(ibasis) = refpot(ibasis)
      kaoezold(1, ibasis) = kaoez(1, ibasis)
      rbasissave(:, ibasis) = rbasis(:, ibasis)
    end do

    ! -----------------------------------------------------------------
    ! read the scoef file
    ! -----------------------------------------------------------------
    open (unit=32452345, file=i25, iostat=ierr)
    write (1337, *) '*', i25, '*'
    write (1337, *) '*', ierr, '*'
    if (ierr/=0) stop '[addvirtual] file not found'
    read (32452345, *) natomimp
    write (1337, *) 'natomimp', natomimp
    allocate (ratomimp(3,natomimp))
    allocate (atomimp(natomimp))
    do iatom = 1, natomimp
      read (32452345, *) ratomimp(:, iatom), atomimp(iatom)
      write (1337, '(A8,I5,A2,3F25.16)') 'IMPATOM ', iatom, ' :', ratomimp(:, iatom)
    end do

    ! -----------------------------------------------------------------
    ! set bulk/surface
    ! -----------------------------------------------------------------
    if (linterface) then
      ndim = 2
      write (1337, '(23X,A)') 'ADDVIRTUAL : surface geometry mode'
    else
      ndim = 3
      write (1337, '(23X,A)') 'ADDVIRTUAL : bulk geometry mode'
    end if


    ! invert bravais vectors
    allocate (bravaisinv(ndim,ndim))
    bravaisinv = bravais(1:ndim, 1:ndim)
    call inverse_d1(bravaisinv(:,:), ndim)


    nrefold = 0
    do ibasis = 1, naez + nemb
      nrefold = max(nrefold, refpotold(ibasis))
    end do
    write (1337, *) 'Number of reference potentials is currently', nrefold

    ! Change basis vectors to cartesian coordinates in order to calculate
    ! distances.
    if (.not. lcartesian) then

      if (linterface) then
        do i = 1, naez + nemb
          do j = 1, ndim
            rbasisold(j, i) = (rbasis(1,i)*bravais(j,1)+rbasis(2,i)*bravais(j,2))
          end do
          rbasisold(3, i) = rbasis(3, i)
        end do
      else
        do i = 1, naez + nemb
          do j = 1, ndim
            rbasisold(j, i) = (rbasis(1,i)*bravais(j,1)+rbasis(2,i)*bravais(j,2)+rbasis(3,i)*bravais(j,3))
          end do
        end do
      end if

    else

      do i = 1, naez + nemb
        rbasisold(:, i) = rbasis(:, i)
      end do

    end if


    ! If the 1st imp. atom in the list is at (0,0,0) then all coordinates are
    ! assumed
    ! relative to the 1st imp atom, otherwise relative to the lattice coords
    ! (absolute coords).
    labscord = .false.
    j = 0
    do j = 1, 3
      if (abs(ratomimp(j,1))>1e-8_dp) labscord = .true.
    end do

    allocate (rclsnew(3,natomimp))
    allocate (rbasisnew1(3,natomimp))
    do i = 1, natomimp
      call dcopy(3, ratomimp(1,i), 1, rclsnew(1,i), 1)
    end do
    if (.not. labscord) then
      iq = atomimp(1)
      do i = 1, natomimp
        call daxpy(3, 1e0_dp, rbasisold(1,iq), 1, rclsnew(1,i), 1)
      end do
    end if
    rmaxclus = 0e0_dp
    do i = 2, natomimp
      diff = 0e0_dp
      do j = 1, 3
        diff = diff + (rclsnew(j,i)-rclsnew(j,1))**2
      end do
      diff = sqrt(diff)
      rmaxclus = max(rmaxclus, diff)
    end do

    nbr(1:3) = 0
    call getclusnxyz(rmaxclus, bravais, ndim, diff, nbr)
    nmax = max(nbr(1), nbr(2), nbr(3))
    nmaxz = nmax
    if (ndim==2) nmaxz = 0
    iq = 0
    do n1 = -nmax, nmax
      do n2 = -nmax, nmax
        do n3 = -nmaxz, nmaxz

          vec1(1:3) = real(n1, kind=dp)*bravais(1:3, 1) + real(n2, kind=dp)*bravais(1:3, 2) + real(n3, kind=dp)*bravais(1:3, 3)

          do i1 = 1, naez
            iq = iq + 1
            diff = 0e0_dp
            vec2(1:3, iq) = vec1(1:3) + rbasisold(1:3, i1)
          end do

        end do
      end do
    end do

    ibasis = 0
    do i = 1, natomimp
      do i1 = 1, iq
        diff = sqrt((rclsnew(1,i)-vec2(1,i1))**2+(rclsnew(2,i)-vec2(2,i1))**2+(rclsnew(3,i)-vec2(3,i1))**2)
        if (diff<=(tol)) go to 100 ! Position is on lattice, do not set as virtual atom
      end do
      call rtobasis(bravais, rclsnew(:,i), rbasisnew, ndim)
      if (linterface) then
        do j = 1, 2
          rbasisnew1(j, i) = rbasisnew(1)*bravaisinv(j, 1) + rbasisnew(2)*bravaisinv(j, 2)
        end do
        rbasisnew1(3, i) = rbasisnew(3)
      else
        rbasisnew1(1:3, i) = rbasisnew(1)*bravaisinv(1:3, 1) + rbasisnew(2)*bravaisinv(1:3, 2) + rbasisnew(3)*bravaisinv(1:3, 3)
      end if
      write (1337, *) 'rnew', rbasisnew1(:, i)
      if (i>1) then
        do i1 = 1, i - 1
          diff = sqrt((rbasisnew1(1,i)-rbasisnew1(1,i1))**2+(rbasisnew1(2,i)-rbasisnew1(2,i1))**2+(rbasisnew1(3,i)-rbasisnew1(3,i1))**2)
          if (diff<=1e-05_dp) go to 100
        end do
      end if
      ibasis = ibasis + 1
      ivir(ibasis) = i
100 end do

    ! IBASIS is the number of virtual atoms
    write (1337, *) 'ibasis', ibasis, (ivir(j), j=1, ibasis)

    if (ibasis+naez>naezd) then
      write (*, *) '[addvirtual] naez increased to ', ibasis
      write (*, *) '[addvirtual] naezd is', naezd
      write (*, *) '[addvirtual] naeznew > naezd please change naezd'
      stop 'addvistual'
    else
      write (1337, *) 'NAEZ will soon be increased to : ', ibasis + naez
    end if


    nvirt = ibasis
    nineq = nineq + nvirt
    ncls = ncls + nvirt
    naeznew = nvirt + naez
    if (naeznew>naezd) then
      write (*, *) '[addvirtual] naez increased to ', naeznew
      write (*, *) '[addvirtual] naezd is', naezd
      write (*, *) '[addvirtual] naeznew > naezd please change naezd'
      stop 'addvirtual'
    end if

    do i = naeznew + 1, naeznew + nemb ! Added +1 : Phivos 25.7.2014
      rbasis(:, i) = rbasissave(:, i-ibasis)
      refpot(i) = refpotold(i-ibasis)
      kaoez(1, i) = kaoezold(1, i-ibasis)
    end do
    do i = naeznew + nemb, naeznew + 1, -1
      rmtrefat(i) = rmtrefat(i-nvirt) ! Shift values of embedded positions in
      ! array
    end do
    do i = naeznew - nemb, naeznew
      rmtrefat(i) = 1.e-20_dp
    end do


    do i = naez + 1, naeznew
      rbasis(:, i) = rbasisnew1(:, ivir(i-naez))
      refpot(i) = nrefold + 1
      noq(i) = 0
      kaoez(1, i) = -1
    end do
    do i = 1, naeznew + nemb
      nref = max(nref, refpot(i))
    end do

    ! -----------------------------------------------------------------
    ! write out stuff
    ! -----------------------------------------------------------------
    write (1337, *) 'addvirtual: List of new basis atoms including virtual atoms'
    do j = 1, naeznew + nemb
      write (1337, *) rbasis(:, j)
    end do
    write (1337, *) '-------------------------------------------------'

    write (1337, *) 'naeznew is now ', naeznew
    write (1337, *) 'setting naez to naeznew'
    naez = naeznew
    write (1337, *) 'updating rbasis array with virtual basis sites'

    do ibasis = 1, naeznew + nemb
      write (1337, *) 'REFPOT', refpot(ibasis)
      write (1337, *) 'NOQ', noq(ibasis)
      write (1337, *) 'KAOEZ', kaoez(1, ibasis)
    end do

    ! stop 'end of ADDVIRTUAL'
    ! deallocate(rbasislist,ratomimp,atomimp)
    deallocate (ratomimp, atomimp)
    deallocate (bravaisinv)
    deallocate (rclsnew)
    deallocate (rbasisnew1)
  end subroutine addviratoms14

  !-------------------------------------------------------------------------------
  !> Summary: check if vec is alreadu in veclist
  !> Author: 
  !> Category: KKRhost, geometry
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Checks if the vector vec is in the vector list veclist in the range of (1,bound)
  !-------------------------------------------------------------------------------
  logical function vec_in_list(vec, veclist, bound)
    use :: mod_datatypes, only: dp
    integer :: bound
    real (kind=dp) :: vec(3)
    real (kind=dp) :: veclist(3, bound)
    integer :: ilist
    real (kind=dp) :: tempvec(3), diff

    vec_in_list = .false.
    do ilist = 1, bound
      tempvec = vec - veclist(:, ilist)
      diff = sqrt(tempvec(1)**2+tempvec(2)**2+tempvec(3)**2)
      if (diff<10e-5_dp) vec_in_list = .true.
    end do
  end function vec_in_list


  !-------------------------------------------------------------------------------
  !> Summary: Converts spacial vector to basis vector
  !> Author: 
  !> Category: KKRhost, geometry
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Converts a spacial vector rpos to a basis vector rbasis such that rbasis = bravais * n with n in [0,1]^ndim
  !-------------------------------------------------------------------------------
  subroutine rtobasis(bravais, rpos, rbasis, ndim)
    use :: mod_datatypes, only: dp
    implicit none
    real (kind=dp), intent (in) :: bravais(3, 3)
    real (kind=dp), intent (in) :: rpos(3)
    integer, intent (in) :: ndim
    real (kind=dp), intent (out) :: rbasis(3)

    real (kind=dp) :: bravais_inv(ndim, ndim)
    real (kind=dp) :: ncoeffreal(ndim)
    integer :: ncoeffint(ndim)
    integer :: idim

    ! --------------------------
    ! first invert the bravais matrix => bravais_inv
    ! --------------------------
    bravais_inv = bravais(1:ndim, 1:ndim)
    call inverse_d1(bravais_inv(:,:), ndim)
    ! --------------------------
    ! then do n = Bravais* rpos
    ! --------------------------
    ncoeffreal = 0.0e0_dp
    do idim = 1, ndim
      ncoeffreal = ncoeffreal + bravais_inv(1:ndim, idim)*rpos(idim)
    end do

    ! ------old method ---------
    ! take the smaller integer value of n => ncoeffint
    ! --------------------------
    !do idim = 1, ndim
    !   ncoeffint(idim) = floor(ncoeffreal(idim))
    !end do

    ! ------new method---------
    ! take the smaller integer value of n + 0.5 => ncoeffint
    ! --------------------------
    do idim = 1, ndim
      ncoeffint(idim) = floor(ncoeffreal(idim)+0.5_dp)
    end do
    ! --------------------------
    ! do rbasis = rpos - bravias * ncoeffint
    ! --------------------------
    rbasis = 0.0e0_dp
    do idim = 1, ndim
      rbasis = rbasis + bravais(:, idim)*ncoeffint(idim)
    end do
    rbasis = rpos - rbasis
  end subroutine rtobasis

  !-------------------------------------------------------------------------------
  !> Summary: Invert matrix
  !> Author: 
  !> Category: KKRhost, undefined
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Computes inverse of a double precision matrix of size n x n using LAPACK calls dgetrf, dgetri 
  !-------------------------------------------------------------------------------
  subroutine inverse_d1(mat, n)
    use :: mod_datatypes, only: dp
    implicit none
    integer, intent (in) :: n
    real (kind=dp), dimension (n, n), intent (inout) :: mat
    real (kind=dp), allocatable :: work(:)
    integer, allocatable :: ipiv(:)
    integer :: info

    allocate (ipiv(n), work(n), stat=info)
    if (info/=0) stop 'error allocating work arrays in inverse_d1 of addviratom'
    call dgetrf(n, n, mat, n, ipiv, info)
    if (info/=0) stop 'inverse_d1: dpotrf failed.'
    call dgetri(n, mat, n, ipiv, work, n, info)
    if (info/=0) stop 'inverse_d1: dpotri failed.'
    deallocate (ipiv, work, stat=info)
    if (info/=0) stop 'error allocating work arrays in inverse_d1 of addviratom'
  end subroutine inverse_d1

end module mod_addvirtual14
