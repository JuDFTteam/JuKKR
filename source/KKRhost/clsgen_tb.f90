!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------
!> Summary: Screening clusters around each atom
!> Author: 
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!> Module providing clsgen_tb subroutine that is used to generate screening 
!> (tight-binding) clusters around the atoms
!-------------------------------------------------------------------------------
module mod_clsgen_tb

  private
  public :: clsgen_tb

contains


  !-------------------------------------------------------------------------------
  !> Summary: Calculates screening cluster of each atom
  !> Author: 
  !> Category: KKRhost, reference-system, geometry
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> This subroutine is used to create the clusters around each atom
  !> (Based on clsgen99.f). Also the reference potential height and radius is
  !> set (vref and rmtref).
  !>
  !> STRATEGY :
  !> Calculate the cluster of each atom by the lattice
  !> parameters avaliable. Sort the atoms in a unique way : big r, big z, big y
  !> compare the positions with the previous clusters to see if there is
  !> a difference. If not keep only previous clusters and make indexing if
  !> a new cluster is found then check dimensions and continue for the new
  !> atom.
  !-------------------------------------------------------------------------------
  subroutine clsgen_tb(naez, nemb, nvirt, rr, rbasis, kaoez, zat, cls, ncls, nacls, atom, ezoa, nlbasis, nrbasis, nleft, nright, zperleft, zperight, tleft, tright, rmtref, &
    rmtrefat, vref, irefpot, nrefpot, rcls, rcut, rcutxy, alat, natyp, nclsd, nrd, naclsd, nrefd, nembd, linterface, nprincd, nprinc, invmod)
    ! ************************************************************************
    ! 
    use :: mod_datatypes, only: dp
    use :: mod_runoptions, only: write_tb_coupling, disable_print_serialnumber
    use :: mod_version_info, only: version_print_header
    use :: mod_dsort, only: dsort
    use :: godfrin, only: t_godfrin ! GODFRIN MdSD
    implicit none
    ! .. arguments
    integer :: naez                !! number of atoms in EZ
    integer :: nemb                !! number of embedding postions
    integer :: ncls                !! number of diff. clusters
    integer :: nlr                 !! =NEMB in decimation, =0 in slab or bulk
    integer :: nvirt               !! Number of virtual atoms
    integer :: nclsd
    integer :: nrd
    integer :: naclsd
    integer :: nrefd
    integer, intent(in) :: nprincd !! user-supplied number of layers in a principal layer
    integer :: nprinc              !! Calculated number of layers in a principal layer
    integer, intent(in) :: invmod  !! MdSD: to check if generalized block tridiagonal inversion is used
    integer :: natyp
    integer :: nembd
    real (kind=dp) :: alat         !! lattice constant A
    real (kind=dp) :: rcut, rcutxy
    real (kind=dp) :: rbasis(3, naez+nembd) !! pos. of basis atoms in EZ
    real (kind=dp) :: rcls(3, naclsd, nclsd) !! real space position of atom in cluster
    real (kind=dp) :: rr(3, 0:nrd) !! set of lattice vectors
    real (kind=dp) :: zat(natyp)   !! nucleus charge
    real (kind=dp) :: rmtref(nrefd)
    real (kind=dp) :: vref(nrefd)

    logical, intent (in) :: linterface


    integer :: cls(naez+nembd)     !! type of cluster around atom
    integer :: kaoez(natyp, naez+nembd) !! type of atom at position in EZ
    integer :: nacls(nclsd)        !! number of atoms in cluster
    integer :: atom(naclsd, naez+nembd) !! index to atom in elem/cell at site in cluster
    integer :: ezoa(naclsd, naez+nembd) !! index to bravais lattice  at site in cluster

    ! .. locals
    ! IX,
    integer :: ilay, n1, ir, isite, jsite, iat1, na, number, maxnumber, pos, ia, in, ib, nb, ii, jatom, icu, ic, iat, i1, icluster, nclsall
    integer :: iatom(naclsd), iezoa(naclsd), isort(naclsd), icouplmat(naez, naez), irep(nclsd) ! representative atom of cluster (inverse of CLS)
    integer :: irefpot(naez+nembd), nrefpot
    real (kind=dp) :: rmtrefat(naez+nembd), rmtref1(naez+nembd)
    real (kind=dp) :: vrefat(naez+nembd), vref1(naez+nembd)


    real (kind=dp) :: r2, epsshl, tol, tol2, distmin, rcls1(3, naclsd), rg(3, naclsd), tmp(3), rsort(naclsd), zat_tmp
    integer :: nlbasis, nrbasis, nleft, nright
    real (kind=dp) :: zperleft(3), zperight(3), tleft(3, nlbasis), tright(3, nrbasis)
    real (kind=dp) :: rcut2, rcutxy2, rxy2, dist

    logical :: lfound

    data epsshl/1.e-4_dp/
    data tol/1.e-7_dp/
    data tol2/1.e-7_dp/


    ! ------------------------------------------------------------------------
    write (1337, *) '>>> CLSGEN_TB: generation of cluster coordinates'
    ! This is generating the clusters which have a distance smaller
    ! than RCUT and RCUTXY in plane .
    ! The cluster atoms are ordered with radius and then z>y>x
    ! The ordering allows an easy comparison of clusters
    ! The principal layer for each layer (atom in unit cell) is
    ! calculated also for each cluster and the maximum number
    ! is returned. Some dimension tests are also done
    write (1337, *) 'RCUT = ', rcut, ' RCUTXY = ', rcutxy
    if (abs(rcutxy-rcut)<1.e-4_dp) then
      write (1337, *) 'Spherical Clusters are created'
      ! LSPHER=.TRUE.
    end if

    ! open clusters file
    open (8, file='clusters', status='unknown')
    call version_print_header(8, disable_print=disable_print_serialnumber)
    write (8, 230) naez
    write (8, 260) alat
    write (8, 240)(zat(kaoez(1,iat)), iat=1, naez-nvirt)
    write (8, 250)(kaoez(1,iat), iat=1, naez)

    rcutxy2 = rcutxy**2
    rcut2 = rcut**2
    nlr = 0
    if (linterface) nlr = nemb
    vrefat(:) = 8.e0_dp            ! Set to 8 Rydbergs
    vref1(:) = 8.e0_dp             ! Set to 8 Rydbergs
    vref(:) = 8.e0_dp              ! Set to 8 Rydbergs

    ! ===============================================================
    ! Check the if the dimension NACLSD is enough before starting to
    ! assign clusters. Write out necessary dimension.
    ! Also find touching MT radius.
    do jatom = 1, naez + nlr       ! loop in all sites incl. left/right host
      ! in decimation case

      maxnumber = 0
      number = 0                   ! counter for atoms in cluster
      distmin = 1.e100_dp          ! Large initial value for RMT**2
      do na = 1, naez              ! loop in all sites in unit cell
        if (kaoez(1,na)/=-1) then  ! Exclude virtual atoms from clusters
          ! (except clust. center)
          do ir = 0, nrd           ! loop in all bravais vectors
            tmp(1:3) = rr(1:3, ir) + rbasis(1:3, na) - rbasis(1:3, jatom)
            rxy2 = tmp(1)**2 + tmp(2)**2
            r2 = tmp(3)**2 + tmp(1)**2 + tmp(2)**2
            if ((rxy2<=rcutxy2) .and. (r2<=rcut2)) then
              number = number + 1
            end if
            if (r2>tol2) distmin = min(distmin, r2)
          end do                   ! IR loop in bravais
        end if                     ! (KAOEZ(1,NA).NE.-1)
      end do                       ! NA loop in NAEZ


      ! In the case of 2 dimensional case loop in the atoms outside.

      if (linterface) then
        ! Somehow messy (lionel messi?)
        do ir = 0, nrd
          do ilay = nleft, 1, -1   ! loop in some layers on left side
            do i1 = nlbasis, 1, -1 ! loop in representative atoms on left side
              tmp(1:3) = rr(1:3, ir) + tleft(1:3, i1) + (ilay-1)*zperleft(1:3) - rbasis(1:3, jatom)
              rxy2 = tmp(1)**2 + tmp(2)**2
              r2 = tmp(3)**2 + rxy2

              if ((rxy2<=rcutxy2) .and. (r2<=rcut2)) then
                number = number + 1
              end if
              if (r2>tol2) distmin = min(distmin, r2)
            end do
          end do

          do ilay = 1, nright
            do i1 = 1, nrbasis
              tmp(1:3) = rr(1:3, ir) + tright(1:3, i1) + (ilay-1)*zperight(1:3) - rbasis(1:3, jatom)
              rxy2 = tmp(1)**2 + tmp(2)**2
              r2 = tmp(3)**2 + rxy2
              if ((rxy2<=rcutxy2) .and. (r2<=rcut2)) then
                number = number + 1
              end if
              if (r2>tol2) distmin = min(distmin, r2)
            end do
          end do


        end do                     ! loop in all bravais lattices
      end if                       ! L2DIM Interface calculation

      distmin = sqrt(distmin)/2.e0_dp ! Touching RMT
      ! Define MT-radius of TB-ref. potentials if not read-in from the input
      ! I.e., if not read in  &
      if (rmtrefat(jatom)<0.e0_dp) rmtrefat(jatom) = int(distmin*alat*99.e0_dp)/100.e0_dp ! round up the decimals
      if (kaoez(1,jatom)==-1) then ! Virtual atom
        rmtrefat(jatom) = 1.e-20_dp
        vrefat(jatom) = 0.e0_dp
      end if


      maxnumber = max(maxnumber, number) ! Find largest cluster size

      write (1337, *) 'clsgen_tb: cluster size of site:', jatom, ':', number
      write (1337, *) 'clsgen_tb: Touching RMT of site:', jatom, ':', distmin


    end do

    if (number>naclsd) then
      write (6, *) '(a) Increase the parameter NACLSD ', 'to a value greater equal ', maxnumber, '.'
      stop 'clsgen_tb: Dimension error (a).'
    end if

    ! ===============================================================

    ! Find different types of ref. potential according to the rmtref and vref
    irefpot(1) = 1
    rmtref1(1) = rmtrefat(1)
    vref1(1) = vrefat(1)
    nrefpot = 1
    do isite = 2, naez + nemb
      lfound = .false.
      do jsite = 1, isite - 1
        if (abs(rmtrefat(isite)-rmtrefat(jsite))+abs(vrefat(isite)-vrefat(jsite))<=tol) then
          irefpot(isite) = irefpot(jsite)
          lfound = .true.
        end if
      end do
      if (.not. lfound) then
        nrefpot = nrefpot + 1
        irefpot(isite) = nrefpot
        ! RMTREFAT goes over all sites, RMTREF1 only over all inequivalent ref.
        ! potentials.
        rmtref1(nrefpot) = rmtrefat(isite)
        vref1(nrefpot) = vrefat(isite)
      end if
    end do

    if (nrefpot>nrefd) then
      write (*, *) 'clsgen_tb: NREFPOT.GT.NREFD:', nrefpot, nrefd
      stop 'clsgen_tb: NREFPOT.GT.NREFD'
    end if
    ! Now that the dimension is known, copy to array RMTREF
    do i1 = 1, nrefpot
      if (rmtref(i1)<0.e0_dp) rmtref(i1) = rmtref1(i1)
    end do
    vref(1:nrefpot) = vref1(1:nrefpot)

    icluster = 1
    rcutxy2 = (rcutxy+epsshl)*(rcutxy+epsshl)
    rcut2 = (rcut+epsshl)*(rcut+epsshl)



    ! ===============================================================

    do jatom = 1, naez + nlr       ! loop in all atoms or layers

      cls(jatom) = 0

      number = 0                   ! counter for sites in cluster
      do na = 1, naez              ! loop in all sites
        if (kaoez(1,na)/=-1) then  ! proceed only if the neighbour is not
          ! virtual atom
          do ir = 0, nrd           ! loop in all bravais vectors
            tmp(1:3) = rr(1:3, ir) + rbasis(1:3, na) - rbasis(1:3, jatom)
            rxy2 = tmp(1)**2 + tmp(2)**2
            r2 = tmp(3)**2 + tmp(1)**2 + tmp(2)**2

            if ((rxy2<=rcutxy2) .and. (r2<=rcut2)) then
              number = number + 1
              atom(number, jatom) = na ! store the atom in elem cell
              ezoa(number, jatom) = ir ! store the bravais vector
              rcls1(1:3, number) = tmp(1:3)
            end if
          end do                   ! IR loop in bravais
        end if
      end do                       ! NA loop in NAEZ


      ! In the case of 2 dimensional case loop in the atoms outside.

      if (linterface) then
        ! Somehow messy (eh? lionel?)
        ! Index ATOM gives the kind of atom

        do ir = 0, nrd
          do ilay = nleft, 1, -1   ! loop in some layers on left side
            do i1 = nlbasis, 1, -1 ! loop in representative atoms on left side
              tmp(1:3) = rr(1:3, ir) + tleft(1:3, i1) + (ilay-1)*zperleft(1:3) - rbasis(1:3, jatom)
              rxy2 = tmp(1)**2 + tmp(2)**2
              r2 = tmp(3)**2 + rxy2

              if ((rxy2<=rcutxy2) .and. (r2<=rcut2)) then
                number = number + 1
                atom(number, jatom) = -naez - i1 ! negative values are used in dlke1
                ezoa(number, jatom) = ir ! ILAY,I1 are negative
                rcls1(1:3, number) = tmp(1:3)
              end if
            end do
          end do


          do ilay = 1, nright
            do i1 = 1, nrbasis
              tmp(1:3) = rr(1:3, ir) + tright(1:3, i1) + (ilay-1)*zperight(1:3) - rbasis(1:3, jatom)
              rxy2 = tmp(1)**2 + tmp(2)**2
              r2 = tmp(3)**2 + rxy2
              if ((rxy2<=rcutxy2) .and. (r2<=rcut2)) then
                number = number + 1
                atom(number, jatom) = -naez - nlbasis - i1
                ezoa(number, jatom) = ir
                rcls1(1:3, number) = tmp(1:3)
              end if
            end do
          end do


        end do                     ! loop in all bravais lattices
      end if                       ! L2DIM Interface calculation

      ! Now the atom JATOM has its cluster.
      ! Sort the atoms of the cluster in increasing order.
      ! First by distance, then by z, then by y, then by x.

      if (number>naclsd) then      ! should not hit here, this was checked
        ! earlier
        write (6, *) '(b) Increase the parameter NACLSD ', 'to a value greater equal ', number, '.'
        stop 'clsgen_tb: Dimension error (b).'
      end if

      do ia = 1, number
        rsort(ia) = sqrt(rcls1(1,ia)**2+rcls1(2,ia)**2+rcls1(3,ia)**2)

        rsort(ia) = 1.e9_dp*rsort(ia) + 1.e6_dp*rcls1(3, ia) + 1.e3_dp*rcls1(2, ia) + 1.e0_dp*rcls1(1, ia)
      end do

      call dsort(rsort, isort, number, pos)
      ! Rearange exchange ia with ib
      ! MAP temporarily to another array
      do ia = 1, number
        rg(1:3, ia) = rcls1(1:3, ia)
        iatom(ia) = atom(ia, jatom)
        iezoa(ia) = ezoa(ia, jatom)
      end do
      ! Now use correct order
      do ia = 1, number
        ib = isort(ia)
        rcls1(1:3, ia) = rg(1:3, ib)
        atom(ia, jatom) = iatom(ib)
        ezoa(ia, jatom) = iezoa(ib)
      end do

      ! Now the clusters have a unique sorting and can be compared with
      ! each other Check if ICLUSTER was found previously

      do icu = 1, icluster - 1
        n1 = nacls(icu)
        iat1 = irep(icu)
        if (clustcomp_tb(rcls,irefpot,atom,iat1,icu,n1,rcls1,number,jatom,naclsd)) cls(jatom) = icu
        ! return true if found before
      end do

      if (cls(jatom)==0) then      ! no equivalent found, add new cluster
        nclsall = icluster         ! incl. embedded atoms of left/right host
        if (jatom<=naez) ncls = icluster ! Excl. embedded atoms of left/right
        ! host
        if (nclsall>nclsd) then
          write (6, *) '(c) Increase the parameter NCLSD ', '  to a value greater equal ', icluster, ' .'
          stop 'clsgen_tb: Dimension error (c).'
        end if
        cls(jatom) = icluster
        nacls(icluster) = number
        irep(icluster) = jatom     ! cluster-class is represented by the
        ! cluster around jatom
        do in = 1, number
          rcls(1:3, in, icluster) = rcls1(1:3, in)
          ! WRITE(6,800) JATOM,ATOM(IN,JATOM),EZOA(IN,JATOM),
          ! &                  (RCLS1(IX,IN),IX=1,3),
          ! &              SQRT(RCLS1(1,IN)**2+RCLS1(2,IN)**2+RCLS1(3,IN)**2)
          ! 800           FORMAT(3I5,4F8.4)
        end do
        icluster = icluster + 1
      end if
      ! ******************************************************
    end do                         ! JATOM = 1,NAEZ + NEMB

    ! Now all clusters of all atoms are found.


    write (1337, *) 'Clusters from clsgen_tb:'
    do jatom = 1, naez
      write (1337, 220) jatom, irefpot(jatom), rmtrefat(jatom), vrefat(jatom), cls(jatom), nacls(cls(jatom))
    end do
    if (linterface) then
      write (1337, *) 'Clusters from clsgen_tb in outer region, left:'
      do ia = 1, nlbasis
        jatom = naez + ia
        write (1337, 220) jatom, irefpot(jatom), rmtrefat(jatom), vrefat(jatom), cls(jatom), nacls(cls(jatom))
      end do

      write (1337, *) 'Clusters from clsgen_tb in outer region, right:'
      do ia = 1, nrbasis
        jatom = naez + nlbasis + ia
        write (1337, 220) jatom, irefpot(jatom), rmtrefat(jatom), vrefat(jatom), cls(jatom), nacls(cls(jatom))
      end do
    end if

    ! Write out clusters in file (already open)
    do jatom = 1, naez
      ic = cls(jatom)
      number = nacls(ic)
      write (8, fmt=150) number
      write (8, fmt=150) jatom, ic
      do i1 = 1, number
        dist = sqrt(rcls(1,i1,ic)**2+rcls(2,i1,ic)**2+rcls(3,i1,ic)**2)
        if (abs(atom(i1,jatom))>natyp) then
          zat_tmp = 0.0_dp
        else
          zat_tmp = zat(abs(atom(i1,jatom)))
        end if
        write (8, 170)(rcls(ii,i1,ic), ii=1, 3), atom(i1, jatom), zat_tmp, dist
      end do

    end do

    ! Write out the coupling matrix
    write (1337, *) 'Coupling matrix:'
    do jatom = 1, naez
      do iat = 1, naez
        icouplmat(jatom, iat) = 0
        do i1 = 1, number
          if (atom(i1,jatom)==iat) then
            icouplmat(jatom, iat) = 1
          end if
        end do
      end do
      write (1337, 280) jatom, (icouplmat(jatom,iat), iat=1, naez)
    end do

    ! ----------------------------------------------------------------------
    if (write_tb_coupling) then
      ! output coupling matrix to file                            ! GODFRIN
      open (file='couplings.dat', unit=123456, status='replace') ! GODFRIN
      write (123456, '("# Couplings between atoms via gref")') ! GODFRIN
      write (123456, '(i8)') naez  ! GODFRIN
      do ia = 1, naez              ! GODFRIN
        write (123456, '(1000i1)') icouplmat(ia, 1:naez) ! GODFRIN
      end do                       ! GODFRIN
      close (123456)               ! GODFRIN
    end if!write_tb_coupling
    ! ----------------------------------------------------------------------

    ! Calculate number of layers in principal layer
    nprinc = 1
    ! MdSD: except if using godfrin
    if (linterface .and. invmod/=3) then
      do jatom = 1, naez           ! loop over rows
        do iat = 1, jatom - 1      ! loop over columns before the diagonal
          if (icouplmat(jatom,iat)==1) nprinc = max(nprinc, jatom-iat)
        end do
        do iat = jatom + 1, naez   ! loop over columns after the diagonal
          if (icouplmat(jatom,iat)==1) nprinc = max(nprinc, iat-jatom)
        end do
      end do
      write (1337, *) 'CLSGEN_TB: Number of layers in a principal layer: NPRINC=', nprinc
    end if

    ! Check if the user-specified blocking is correct
    if (mod(naez,nprincd) == 0) then
      ! MdSD: generalized blocking
      if (invmod == 3) then
        nb = t_godfrin%nb
        do ib=1,nb-1
!         diagonal blocks
          icouplmat(1+sum(t_godfrin%bdims(1:ib-1)):sum(t_godfrin%bdims(1:ib)),1+sum(t_godfrin%bdims(1:ib-1)):sum(t_godfrin%bdims(1:ib))) = 0
!         subdiagonal blocks
          icouplmat(1+sum(t_godfrin%bdims(1:ib)):sum(t_godfrin%bdims(1:ib+1)),1+sum(t_godfrin%bdims(1:ib-1)):sum(t_godfrin%bdims(1:ib))) = 0
!         superdiagonal blocks
          icouplmat(1+sum(t_godfrin%bdims(1:ib-1)):sum(t_godfrin%bdims(1:ib)),1+sum(t_godfrin%bdims(1:ib)):sum(t_godfrin%bdims(1:ib+1))) = 0
        end do
!       last diagonal blocks
        icouplmat(1+sum(t_godfrin%bdims(1:nb-1)):sum(t_godfrin%bdims(1:nb)),1+sum(t_godfrin%bdims(1:nb-1)):sum(t_godfrin%bdims(1:nb))) = 0
!       additional blocks from periodic boundary conditions
        if (t_godfrin%lper) then
          icouplmat(1+sum(t_godfrin%bdims(1:0)):sum(t_godfrin%bdims(1:1)),1+sum(t_godfrin%bdims(1:nb-1)):sum(t_godfrin%bdims(1:nb))) = 0
          icouplmat(1+sum(t_godfrin%bdims(1:nb-1)):sum(t_godfrin%bdims(1:nb)),1+sum(t_godfrin%bdims(1:0)):sum(t_godfrin%bdims(1:1))) = 0
        end if
      ! MdSD: standard blocking
      else
        nb = naez/nprincd
        do ib=1,nb-1
!         diagonal blocks
          icouplmat(1+(ib-1)*nprincd:ib*nprincd,1+(ib-1)*nprincd:ib*nprincd) = 0
!         subdiagonal blocks
          icouplmat(1+ib*nprincd:(ib+1)*nprincd,1+(ib-1)*nprincd:ib*nprincd) = 0
!         superdiagonal blocks
          icouplmat(1+(ib-1)*nprincd:ib*nprincd,1+ib*nprincd:(ib+1)*nprincd) = 0
        end do
!       last diagonal blocks
        icouplmat(1+(nb-1)*nprincd:nb*nprincd,1+(nb-1)*nprincd:nb*nprincd) = 0
      end if
!     mistakes were made?
      if (any(icouplmat /= 0) .and. invmod /= 0) then
        write (1337, *) 'CLSGEN_TB: WARNING --- User-supplied NPRINCD=', nprincd, ' is INCORRECT!'
        ! Write out the coupling matrix
        write (1337, *) 'Truncated coupling matrix:'
        do jatom = 1, naez
          write (1337, 280) jatom, (icouplmat(jatom,iat), iat=1, naez)
        end do
      else
        nprinc = nprincd
      end if
    else
    ! MdSD: error message here?
    end if

    ! close clusters file
    close (8)

    ! ------------------------------------------------------------------------
    write (1337, *) ' Sub clsgen_tb  exiting <<<<<<<<<<<<<'
    ! ------------------------------------------------------------------------

100 format (' cluster around atom     ', 10i4, /, ' number atoms in cluster ', 10i4)
110 format (i4, 2i5, 3f8.2, i6, 4f7.2)
120 format (' cocls : naez =', i3, ' (x,y,z)= ', 3f10.4)
130 format (12x, i6, 3f10.4)
140 format ('  Nr  naez kaoez     x       y       z', '  ezoa  RR(1)  RR(2)  RR(3)      R')
150 format (3i8)
160 format ('Cu  ', 3d24.16, 2i8, f18.12)
170 format (3e27.19, i5, f5.1, e17.9)
180 format (3f12.7, '  scaling factor')
190 format (i4, 3f12.7, '  center', /, (i4,3f12.7))
200 format (i4, 3f12.7, '  center of gravity')
210 format ('contains ', i4, '  atoms.')
220 format ('CLSGEN_TB: Atom', i5, ' Refpot', i5, ' Rmtref', f10.7, ' Vref', f10.7, ' TB-cluster', i5, ' Sites', i5)
230 format (i4)
240 format (('# ZAT   ',20f6.2))
250 format (('# KAOEZ ',20i4))
260 format (f12.7, 6x, 'ALAT')
270 format ('> cluster ', i4, ' at atom ', i4, ' of type ', i4, '.')
280 format (i4, 1x, 500i1)
290 format (i6, 3f15.6)
300 format ('The Number of layers is each Principal Layer = ', i5)


    return
  end subroutine clsgen_tb


  !-------------------------------------------------------------------------------
  !> Summary: Compares two cluster for equivalence
  !> Author: 
  !> Category: KKRhost, geometry, reference-system
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> This function returns true if cluster ic1 is equal to new cluster
  !> RCLS        coordinates of all (already found) clusters
  !> IC1         First cluster index
  !> N1          Number of atoms in IC1 cluster
  !> rcls1       coordinates of new cluster
  !> n2          number of atoms in ic2 cluster
  !> ATOM:       atom-type at a certain position in a cluster
  !> IAT1, IAT2: central atoms of 1st,2nd cluster
  !> IREFPOT:    Type of reference potential
  !-------------------------------------------------------------------------------
  logical function clustcomp_tb(rcls, irefpot, atom, iat1, ic1, n1, rcls1, n2, iat2, naclsd)
    ! 
    use :: mod_datatypes, only: dp
    implicit none
    integer :: naclsd
    real (kind=dp) :: rcls(3, naclsd, *), rcls1(3, naclsd)
    integer :: atom(naclsd, *), irefpot(*)
    integer :: ic1, n1, n2, iat1, iat2, n, i
    real (kind=dp) :: rd, tol
    logical :: lreflog

    tol = 1.e-5_dp
    clustcomp_tb = .false.
    lreflog = .true.
    if (n1==n2) then
      rd = 0.e0_dp
      do n = 1, n1
        ! compare ref-potential types
        if (irefpot(abs(atom(n,iat1)))/=irefpot(abs(atom(n,iat2)))) lreflog = .false.
        do i = 1, 3
          rd = rd + abs(rcls(i,n,ic1)-rcls1(i,n)) ! compare coordinares
        end do
      end do
      if (abs(rd)<tol .and. lreflog) clustcomp_tb = .true.
    end if
  end function clustcomp_tb


end module mod_clsgen_tb
