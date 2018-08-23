module mod_clsgen99

contains

! -------------------------------------------------------------------------------
! SUBROUTINE: CLSGEN99
! > @brief This subroutine is used to create the clusters around each atom
! > where repulsive potentials will be positioned.
! > @details Calculate the cluster of each atom by the lattice parameters
! avaliable.
! > Sort the atoms in a unique way : big r, big z, big y
! > compare the positions with the previous clusters to see if there is
! > a difference. If not keep only previous clusters and make indexing if
! > a new cluster is found then check dimensions and continue for the new
! > atom.
! > @note Small bug in assigning clusters removed 29/04/2003 v.popescu
! > IATCLUS(NCLSD) is pointing to the first atomic site associated with a
! tb-cluster
! > @note JC: This routine seems to be deprecated and not used in the actual
! code
! -------------------------------------------------------------------------------
subroutine clsgen99(naez, rr, nr, rbasis, kaoez, z, cls, nacls, refpot, atom, &
  ezoa, nlbasis, nrbasis, nleft, nright, zperleft, zperight, tleft, tright, &
  rcls, rcut, rcutxy, l2dim, alat, natyp, nemb, nprincd, naclsd, nclsd)

  use :: mod_version_info
  use :: mod_datatypes

  implicit none
  ! .. Arguments ..
  integer, intent (in) :: nr       ! < Number of real space vectors rr
  integer, intent (in) :: nemb     ! < Number of 'embedding' positions
  integer, intent (in) :: naez     ! < Number of atoms in unit cell
  integer, intent (in) :: natyp    ! < Number of kinds of atoms in unit cell
  integer, intent (in) :: nclsd    ! < Maximum number of different TB-clusters
  integer, intent (in) :: nleft    ! < Number of repeated basis for left host
                                   ! to get converged electrostatic potentials
  integer, intent (in) :: nright   ! < Number of repeated basis for right host
                                   ! to get converged electrostatic potentials
  integer, intent (in) :: naclsd   ! < Maximum number of atoms in a TB-cluster
  integer, intent (in) :: nprincd  ! < Number of principle layers, set to a
                                   ! number >= NRPINC in output of main0
  integer, intent (in) :: nlbasis  ! < Number of basis layers of left host
                                   ! (repeated units)
  integer, intent (in) :: nrbasis  ! < Number of basis layers of right host
                                   ! (repeated units)
  real (kind=dp), intent (in) :: alat ! < Lattice constant in a.u.
  real (kind=dp), intent (in) :: rcut ! < Parameter for the screening cluster
                                      ! along the z-direction
  real (kind=dp), intent (in) :: rcutxy ! < Parameter for the screening
                                        ! cluster along the x-y plane
  integer, dimension (naez+nemb), intent (in) :: refpot ! < Ref. pot. card  at
                                                        ! position
  real (kind=dp), dimension (natyp), intent (in) :: zat ! < Nuclear charge
  real (kind=dp), dimension (3, 0:nr), intent (in) :: rr ! < Set of real space
                                                         ! vectors (in a.u.)
  real (kind=dp), dimension (3, naez+nemb), intent (in) :: rbasis ! < Position
                                                                  ! of atoms
                                                                  ! in the
                                                                  ! unit cell
                                                                  ! in units
                                                                  ! of bravais
                                                                  ! vectors
  real (kind=dp), dimension (3, naclsd, nclsd), intent (in) :: rcls ! < Real
                                                                    ! space
                                                                    ! position
                                                                    ! of atom
                                                                    ! in
                                                                    ! cluster

  integer, dimension (naez+nemb), intent (inout) :: cls ! < Cluster around
                                                        ! atomic sites
  integer, dimension (nclsd), intent (inout) :: nacls ! < Number of atoms in
                                                      ! cluster
  integer, dimension (naclsd, naez+nemb), intent (in) :: atom ! < Atom at site
                                                              ! in cluster
  integer, dimension (naclsd, naez+nemb), intent (in) :: ezoa ! < EZ of atom
                                                              ! at site in
                                                              ! cluster
  integer, dimension (natyp, naez+nemb), intent (in) :: kaoez ! < Kind of atom
                                                              ! at site in
                                                              ! elem. cell
  ! .. Local variables

  integer :: i, n1, inum, isum, na, number, n, nprin, itest1, itest
  integer :: pos, ia, in, ib, ii, jatom, icu, ic, iat, i0, i1, icluster
  integer, dimension (naclsd) :: iatom
  integer, dimension (naclsd) :: iezoa
  integer, dimension (naclsd) :: isort
  integer, dimension (nclsd) :: iatcls
  integer, dimension (naez, naez) :: icouplmat

  real (kind=dp) :: r, r2, epsshl
  real (kind=dp) :: rcut2, rcutxy2, rxy2
  real (kind=dp), dimension (3) :: tmp
  real (kind=dp), dimension (naclsd) :: rsort
  real (kind=dp), dimension (3) :: zperleft
  real (kind=dp), dimension (3) :: zperight
  real (kind=dp), dimension (3, naclsd) :: rg
  real (kind=dp), dimension (3, naclsd) :: rcls1
  real (kind=dp), dimension (3, nemb+1) :: tleft
  real (kind=dp), dimension (3, nemb+1) :: tright

  logical :: l2dim, clustcomp

  logical :: test, lspher
  external :: dsort, clustcomp
  intrinsic :: min, sqrt

  data epsshl/1.0d-4/
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This is generating the clusters which have a distance smaller
  ! than RCUT and RCUTXY in plane .
  ! The cluster atoms are ordered with radious and then z>y>x
  ! The ordering allows an easy comparison of clusters
  ! The principal layer for each layer (atom in unit cell) is
  ! calculated also for each cluster and the maximum number
  ! is returned. Some dimension tests are also done
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! ---------------------------------------------------------------------------
  ! OUTPUT
  ! ---------------------------------------------------------------------------
  write (1337, '(79(1H=))')
  write (1337, '(16X,A)') 'CLSGEN99: generation of TB-clusters coordinates'
  write (1337, '(79(1H=))')
  write (1337, *)
  ! ---------------------------------------------------------------------------
  ! OUTPUT
  ! ---------------------------------------------------------------------------
  lspher = .false.
  write (1337, *) 'RCUT = ', rcut, ' RCUTXY = ', rcutxy
  if (abs(rcutxy-rcut)<1.d-4) then
    write (1337, *) 'Spherical Clusters are created'
    lspher = .true.
  end if
  ! ----------------------------------------------------------------------------
  if (test('clusters')) then
    open (8, file='clusters', status='UNKNOWN')
    call version_print_header(8)
    write (8, 130) naez
    write (8, 160) alat
    write (8, 140)(z(kaoez(1,i)), i=1, naez)
    write (8, 150)(kaoez(1,i), i=1, naez)
  end if
  ! ----------------------------------------------------------------------------
  icluster = 1
  do n = 1, nclsd
    iatcls(n) = 0
    nacls(n) = 0
  end do
  call rinit(3*naclsd*nclsd, rcls)

  rcutxy2 = (rcutxy+epsshl)*(rcutxy+epsshl)
  rcut2 = (rcut+epsshl)*(rcut+epsshl)

  do jatom = 1, naez               ! loop in all atoms or layers
    cls(jatom) = 0
    number = 0                     ! counter for atoms in cluster
    do na = 1, naez                ! loop in all atoms
      do n = 0, nr                 ! loop in all bravais vectors
        do i = 1, 3
          tmp(i) = rr(i, n) + rbasis(i, na) - rbasis(i, jatom)
        end do
        rxy2 = tmp(1)**2 + tmp(2)**2
        r2 = tmp(3)**2
        if (lspher) r2 = r2 + rxy2

        if ((rxy2<=rcutxy2) .and. (r2<=rcut2)) then
          number = number + 1
          if (number>naclsd) then
            write (6, *) ' ERROR: Dimension NACLSD in inc.cls too small', &
              number, naclsd
            stop '   < CLSGEN99 >'
          end if

          atom(number, jatom) = na ! store the atom in elem cell
          ezoa(number, jatom) = n  ! store the bravais vector
          do i = 1, 3
            rcls1(i, number) = tmp(i)
          end do
        end if
      end do                       ! N loop in bravais
    end do                         ! NA loop in NAEZ

    ! -------------------------------------------------------------------------
    ! In the case of 2 dimensional case loop in the atoms
    ! outside.
    ! -------------------------------------------------------------------------
    if (l2dim) then
      ! ----------------------------------------------------------------------
      ! Somehow meshy
      ! ATOM gives the kind of atom
      ! ----------------------------------------------------------------------
      do n = 0, nr
        do i = nleft, 1, -1        ! loop in some layers on left side
          do i1 = nlbasis, 1, -1   ! loop in representative atoms on left side
            do i0 = 1, 3
              tmp(i0) = rr(i0, n) + tleft(i0, i1) + (i-1)*zperleft(i0) - &
                rbasis(i0, jatom)
            end do
            rxy2 = tmp(1)**2 + tmp(2)**2
            r2 = tmp(3)**2
            if (lspher) r2 = r2 + rxy2

            if ((rxy2<=rcutxy2) .and. (r2<=rcut2)) then

              number = number + 1
              if (number>naclsd) then
                write (6, *) 'ERROR: Dimension NACLSD in inc.cls too small', &
                  number,

end module mod_clsgen99
