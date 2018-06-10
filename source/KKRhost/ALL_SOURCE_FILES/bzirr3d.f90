subroutine bzirr3d(nkp, nkxyz, kpoibz, kp, recbv, bravais, wtkp, volbz, &
  rsymat, nsymat, isymindex, symunitary, irr, krel, iprint)
!===========================================================================
!info
!info   find irreducible BZ and create mesh in it.
!info           original version Arthur Ernst, Daresbury, 27/03/98
!                                fcc lattice bug corrected
!                                modified 26/5/99, 2/06/99
! Modified on 20.01.2000 To use the 2d inversion!
! Modified Dec 2001 - Apr 2002 to deal with relativistic case, bugs
!          removed (full BZ integration gives now same results as
!          symetrised calculation)   HE/VP, Munich

!===========================================================================
!  Input:
!         nkxyz : original k-mesh net in the 3 directions of the reciprocal
!                 lattice vectors (not xyz directions).
!         recbv : reciprocal unit cell
!                 (normalized: recbv(i,j)=G(i,j)*alat/(2*pi), i=x,y,z j=1,2,3)
!                 if G is the "normal" rec. lattice vector.
!      bravais  : direct unit cell in a.u.
!      rsymat   : symmetry rotation matrices of real lattice
!      nsymat   : number of rotation matrices of real lattice
!         irr   : if .true. then irreducible BZ else take all BZ

!     symunitary: in case of a FM REL calculation, indicates a unitary/
!                 antiunitary symmetry operation

!  Output:
!         nkp   : number of points in irreducible BZ
!         kp    : k-point mesh
!         wtkp  : weights for k-points
!        volbz  : volume of the BZ
!  Inside:
!         ibk   : Flag showing if the mesh point jx,jy,jz has already been
!                 taken. Could also be a logical variable.
!==========================================================================
  implicit none
  integer :: maxk1, maxk2, maxk3, nsymaxd
  parameter (maxk1=501, maxk2=501, maxk3=100, nsymaxd=48)
!PARAMETER (MAXK1=350,MAXK2=350,MAXK3=70,NSYMAXD=48)
! i/o
! made local ibk array allocatable to not take too much memory
!integer nkp,nkxyz(3),ibk(0:MAXK1,0:MAXK2,0:MAXK3),kpoibz,nsymat
  integer :: nkp, nkxyz(3), kpoibz, nsymat
  integer :: isymindex(*), nkxyz1(3), krel
  double precision :: kp(3, *), wtkp(*), volbz, cf(3)
  double precision :: recbv(3, 3), bravais(3, 3), rsymat(64, 3, 3)
  logical :: irr
! local
  integer, allocatable :: ibk(:, :, :)
  real *8 :: bginv(3, 3), bgmat(3, 3), bgp(3, 3), bv(3)
  integer :: nbgp(3, 3, 48), isym
  integer :: iktab(3), ind2(3), iprint
  integer :: i, j, jx, jy, jz, nsym, iws, iwt, is, ja, jb, jc, ix, iy, iz, k, &
    nk, n
  integer :: ndim
  double precision :: u(3, 3, 48), gq(3, 3), v1
  logical :: lsurf
  logical :: symunitary(nsymaxd)
  double precision :: msign

  double precision :: ddet33, ddot
  external :: ddet33, ddot

  allocate (ibk(0:maxk1,0:maxk2,0:maxk3))


!heck if we are in surface mode

  lsurf = .false.
  if (bravais(1,3)==0.d0 .and. bravais(2,3)==0.d0 .and. bravais(3,3)==0.d0) &
    lsurf = .true.

  ndim = 3
  if (lsurf) then
    nkxyz(3) = 1
    ndim = 2
  end if

  if (nkxyz(1)>maxk1 .or. nkxyz(2)>maxk2 .or. nkxyz(3)>maxk3) then
    write (6, *) 'BZIRR3D : Increase MAXK ', (nkxyz(i), i=1, 3), maxk1, maxk2, &
      maxk3
    stop
  end if
  do i = 1, 3
    nkxyz1(i) = nkxyz(i)
  end do
!-------------------

! Create small unit cell for the integration in the reciprocal space (gq),

!                                          steps along the basis vectors
  do i = 1, 3
    do j = 1, 3
      gq(i, j) = recbv(i, j)/dble(nkxyz1(j))
    end do
  end do
!=======================================================================
  if (irr) then

    nsym = nsymat

!-----------------------------------------------------------------------

    if (nsym>48) stop 'DIM problem in bzirr3d.'

    if ((lsurf) .and. (nsym>12)) stop 'DIM problem in bzirr3d, surf mode.'

!----------------------------------------------------------------------

    do n = 1, nsym
      isym = isymindex(n)
      msign = 1.d0
      if ((krel==1) .and. (.not. symunitary(n))) msign = -1.d0
      do i = 1, 3
        do j = 1, 3
          u(i, j, n) = msign*rsymat(isym, i, j)
        end do
      end do
    end do
!========================================================================
  else
!========================================================================
    nsym = 1
    u(1, 1, 1) = 1.d0
    u(1, 2, 1) = 0.d0
    u(1, 3, 1) = 0.d0
    u(2, 1, 1) = 0.d0
    u(2, 2, 1) = 1.d0
    u(2, 3, 1) = 0.d0
    u(3, 1, 1) = 0.d0
    u(3, 2, 1) = 0.d0
    u(3, 3, 1) = 1.d0
  end if
!==========================================================================

  do j = 1, 3
    do i = 1, 3
      bgmat(i, j) = ddot(3, gq(1,i), 1, gq(1,j), 1)
    end do
  end do

  call rinvgj(bginv, bgmat, 3, ndim)

!-----------------------------------------------------------------------
!          rotate the 3 step vectors   GQ  --  it should be possible
!          to express the rotated vectors  B(j) in terms of the old ones
!     B(i) = SUM(j) GQ(j) * n(j,i)  with integer coefficients n(j,i)

!==========================================================================
  do is = 1, nsym
    if (iprint>2) write (1337, 100) is

    do i = 1, 3
      call dgemv('N', 3, 3, 1d0, u(1,1,is), 3, gq(1,i), 1, 0d0, bgp(1,i), 1)

      do j = 1, 3
        bv(j) = ddot(3, gq(1,j), 1, bgp(1,i), 1)
      end do
      call dgemv('N', 3, 3, 1d0, bginv, 3, bv, 1, 0d0, cf, 1)

      do j = 1, ndim
        if (abs(nint(cf(j))-cf(j))>1d-8) write (1337, 120) i, j, cf(j)
        nbgp(j, i, is) = nint(cf(j))
      end do
      if (iprint>2) write (1337, 110) i, (bgp(j,i), j=1, 3), bv, cf, &
        (nbgp(j,i,is), j=1, 3)

    end do

  end do
!========================================================================

  nk = nkxyz1(1)*nkxyz1(2)*nkxyz1(3)

  do i = 0, nkxyz1(1)
    do j = 0, nkxyz1(2)
      do k = 0, nkxyz1(3)
        ibk(i, j, k) = 0
      end do
    end do
  end do

  ix = nkxyz1(1)
  iy = nkxyz1(2)
  iz = nkxyz1(3)
  nkp = 0
  iws = 0

  do jx = 0, ix - 1
    iktab(1) = jx
    do jy = 0, iy - 1
      iktab(2) = jy
      do jz = 0, iz - 1
        iktab(3) = jz

        if (ibk(jx,jy,jz)==0) then
          nkp = nkp + 1
          iwt = 0

!========================================================================
          do isym = 1, nsym

!               rotate k-vector  NKP  and transform into parallelepiped
!          ROT k = SUM(i) m(i) ROT gq(i)
!                = SUM(i) m(i) SUM(j) n(i,j) gq(j)
!                = SUM(j) [SUM(i) m(i) n(i,j)] gq(j)

            do j = 1, 3
              is = 0
              do i = 1, 3
                is = is + iktab(i)*nbgp(j, i, isym)
              end do
              is = mod(is, nkxyz1(j))
              if (is<0) is = is + nkxyz1(j)
              ind2(j) = is
            end do

            ja = ind2(1)
            jb = ind2(2)
            jc = ind2(3)

            if (ibk(ja,jb,jc)==0) then
              ibk(ja, jb, jc) = nkp
              iwt = iwt + 1
            end if
! Mesh point in the unit cell found.
          end do
!========================================================================

          if (nkp<=kpoibz) then
            do i = 1, 3
              kp(i, nkp) = 0d0
              do j = 1, 3
                kp(i, nkp) = kp(i, nkp) + gq(i, j)*dble(iktab(j))
              end do
            end do
            wtkp(nkp) = dble(iwt)/dble(nk)
          else
            write (6, *) ' No. of k-points ', nkp, ' > KPOIBZ '
            stop ' <BZIRR3D> increase parameter KPOIBZ'
          end if

        end if
        iws = iws + iwt

      end do
    end do
  end do

  volbz = 0.d0

  if (lsurf) then
    v1 = dabs(recbv(1,1)*recbv(2,2)-recbv(1,2)*recbv(2,1))
  else
    v1 = ddet33(recbv)
  end if

  do i = 1, nkp
    wtkp(i) = wtkp(i)*v1/dble(nsym)
    volbz = volbz + wtkp(i)*dble(nsym)
  end do

  deallocate (ibk)

  return

100 format (5x, 'rotated GQ  for IROT=', i3)
110 format (5x, i3, 3f7.3, 2x, 3f7.3, 2x, 3f7.3, 2x, 3i3)
120 format (5x, 2i3, 3f7.3)
end subroutine
