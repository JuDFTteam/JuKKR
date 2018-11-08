!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_bzkmesh

  private
  public :: bzkmesh

contains

  !-------------------------------------------------------------------------------
  !> Summary: Find different k-meshes
  !> Author: 
  !> Category: KKRhost, k-points
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !> 
  !-------------------------------------------------------------------------------
  subroutine bzkmesh(nbxin, nbyin, nbzin, maxmesh, lirr, bravais, recbv, nsymat, rsymat, isymindex, symunitary, ielast, ez, kmesh, iprint, krel, kpoibz, maxmshd)
use :: mod_runoptions, only: print_kmesh, set_kmesh_large, set_kmesh_small, write_kpts_file, write_rhoq_input --manopt-- 
    use :: mod_types, only: t_inc
    use :: mod_wunfiles, only: t_params
    use :: mod_rhoqtools, only: rhoq_write_kmesh
    use :: mod_datatypes, only: dp
    use :: mod_bzirr3d, only: bzirr3d
    implicit none
    real (kind=dp), parameter :: eps = 1.0d-12
    ! ..
    ! .. Scalar Arguments ..
    integer :: maxmesh, nbxin, nbyin, nbzin, nsymat, iprint, krel, kpoibz, ielast, maxmshd
    logical :: lirr
    ! ..
    ! .. Array Arguments ..
    real (kind=dp) :: bravais(3, 3), recbv(3, 3)
    real (kind=dp) :: rsymat(64, 3, 3)
    integer :: isymindex(*), kmesh(*)
    complex (kind=dp) :: ez(*)
    ! .. unitary/antiunitary symmetry flag
    logical :: symunitary(*)
    ! ..
    ! .. Local Scalars ..
    real (kind=dp) :: volbz
    integer :: i, id, ks, l, n, nbx, nby, nbz, nofks
    ! ..
    ! .. Local Arrays ..
    real (kind=dp) :: bzkp(3, kpoibz), volcub(kpoibz)
    integer :: nxyz(3)
    ! ..
    ! .. External Functions ..
    logical :: test
    external :: test
    ! ---------------------------------------------------------------------

    ! --> set number of different K-meshes

    maxmesh = 1
    if (set_kmesh_large) then
      do i = 1, ielast
        kmesh(i) = 1
      end do
    else
      do i = 1, ielast
        if (abs(aimag(ez(ielast)))>eps) then
          n = int(1.001d0+log(aimag(ez(i))/aimag(ez(ielast)))/log(2.0d0))
        else
          n = 1
        end if
        kmesh(i) = n
        maxmesh = max(maxmesh, n)
        if (kmesh(i)<1) kmesh(i) = 1
      end do
      kmesh(1) = maxmesh
    end if

    if (set_kmesh_small) then
      do i = 1, ielast
        kmesh(i) = maxmesh
      end do
    end if

    if (maxmesh>maxmshd) then
      write (6, fmt='(5X,A,I2)') 'Dimension ERROR: Please increase MAXMSHD to ', maxmesh
      write (6, fmt='(22X,A,/)') 'in the programs < main0 > and < main1b >'
      stop '        < BZKMESH >'
    end if
    ! ---------------------------------------------------------------------
    nbx = nbxin
    nby = nbyin
    nbz = nbzin

    if (write_kpts_file) open (52, file='kpoints', form='formatted')

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    write (1337, 110)
    write (1337, 120) maxmesh, nsymat
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

    ! save maxmesh and allocate kmesh for later use in t_inc and t_params
    t_inc%nkmesh = maxmesh
    t_params%kpoibz = kpoibz
    t_params%maxmesh = maxmesh
    allocate (t_inc%kmesh(maxmesh))
    allocate (t_params%bzkp(3,kpoibz,maxmesh), t_params%volcub(kpoibz,maxmesh), t_params%volbz(maxmesh), t_params%nofks(maxmesh))
    ! needed for wavefunction saving
    allocate (t_inc%kmesh_ie(ielast))
    t_inc%kmesh_ie = kmesh(1:ielast)
    ! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    do l = 1, maxmesh
      if (l>1) then
        nbx = int(nbx/1.4)
        nby = int(nby/1.4)
        nbz = int(nbz/1.4)
      end if
      if (nbx<1) nbx = 1
      if (nby<1) nby = 1
      if (nbz<1) nbz = 1
      nxyz(1) = nbx
      nxyz(2) = nby
      nxyz(3) = nbz

      call bzirr3d(nofks, nxyz, kpoibz, bzkp, recbv, bravais, volcub, volbz, rsymat, nsymat, isymindex, symunitary, lirr, krel, iprint)

      ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      write (1337, 130) l, nofks, (nxyz(i), i=1, 3), volbz
      if (l==maxmesh) write (1337, 140)
      ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

      if (write_kpts_file) write (52, fmt='(I8,F15.10,/,(3F12.8,D20.10))') nofks, volbz, ((bzkp(id,i),id=1,3), volcub(i), i=1, nofks)
      if (write_rhoq_input .and. (l==1)) then
        call rhoq_write_kmesh(nofks, nxyz, volbz, bzkp, volcub, recbv, bravais)
      end if

      t_params%nofks(l) = nofks
      t_params%volbz(l) = volbz
      do i = 1, nofks
        do id = 1, 3
          t_params%bzkp(id, i, l) = bzkp(id, i)
        end do
        t_params%volcub(i, l) = volcub(i)
      end do
      ! save nofks for this mesh in t_inc
      t_inc%kmesh(l) = nofks
      ! ---------------------------------------------------------------------

      ! -->  output of k-mesh

      if (print_kmesh) then
        do ks = 1, nofks
          write (1337, fmt=100)(bzkp(i,ks), i=1, 3), volcub(ks)
        end do
      end if
      ! ---------------------------------------------------------------------
    end do
    ! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    if (write_kpts_file) close (52)
    ! CLOSE (52)
100 format (3f12.5, f15.8)
110 format (5x, '< BZKMESH > : creating k-mesh,', ' write to file kpoints', /)
120 format (8x, 'number of different k-meshes :', i2, /, 8x, 'the direct lattice', i3, ' symmetries will be used', /, /, 8x, 35('-'), /, 8x, 'k-mesh NofKs N kx N ky N kz vol BZ', &
      /, 8x, 35('-'))
130 format (8x, 2i6, 3i5, f8.4)
140 format (8x, 35('-'), /)
  end subroutine bzkmesh

end module mod_bzkmesh
