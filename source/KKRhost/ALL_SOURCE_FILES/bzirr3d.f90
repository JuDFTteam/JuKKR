    Subroutine bzirr3d(nkp, nkxyz, kpoibz, kp, recbv, bravais, wtkp, volbz, &
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
      Use mod_datatypes, Only: dp
      Implicit None
      Integer :: maxk1, maxk2, maxk3, nsymaxd
      Parameter (maxk1=501, maxk2=501, maxk3=100, nsymaxd=48)
!PARAMETER (MAXK1=350,MAXK2=350,MAXK3=70,NSYMAXD=48)
! i/o
! made local ibk array allocatable to not take too much memory
!integer nkp,nkxyz(3),ibk(0:MAXK1,0:MAXK2,0:MAXK3),kpoibz,nsymat
      Integer :: nkp, nkxyz(3), kpoibz, nsymat
      Integer :: isymindex(*), nkxyz1(3), krel
      Real (Kind=dp) :: kp(3, *), wtkp(*), volbz, cf(3)
      Real (Kind=dp) :: recbv(3, 3), bravais(3, 3), rsymat(64, 3, 3)
      Logical :: irr
! local
      Integer, Allocatable :: ibk(:, :, :)
      Real (Kind=dp) :: bginv(3, 3), bgmat(3, 3), bgp(3, 3), bv(3)
      Integer :: nbgp(3, 3, 48), isym
      Integer :: iktab(3), ind2(3), iprint
      Integer :: i, j, jx, jy, jz, nsym, iws, iwt, is, ja, jb, jc, ix, iy, iz, &
        k, nk, n
      Integer :: ndim
      Real (Kind=dp) :: u(3, 3, 48), gq(3, 3), v1
      Logical :: lsurf
      Logical :: symunitary(nsymaxd)
      Real (Kind=dp) :: msign

      Real (Kind=dp) :: ddet33, ddot
      External :: ddet33, ddot

      Allocate (ibk(0:maxk1,0:maxk2,0:maxk3))


!heck if we are in surface mode

      lsurf = .False.
      If (bravais(1,3)==0.E0_dp .And. bravais(2,3)==0.E0_dp .And. &
        bravais(3,3)==0.E0_dp) lsurf = .True.

      ndim = 3
      If (lsurf) Then
        nkxyz(3) = 1
        ndim = 2
      End If

      If (nkxyz(1)>maxk1 .Or. nkxyz(2)>maxk2 .Or. nkxyz(3)>maxk3) Then
        Write (6, *) 'BZIRR3D : Increase MAXK ', (nkxyz(i), i=1, 3), maxk1, &
          maxk2, maxk3
        Stop
      End If
      Do i = 1, 3
        nkxyz1(i) = nkxyz(i)
      End Do
!-------------------

! Create small unit cell for the integration in the reciprocal space (gq),

!                                          steps along the basis vectors
      Do i = 1, 3
        Do j = 1, 3
          gq(i, j) = recbv(i, j)/real(nkxyz1(j), kind=dp)
        End Do
      End Do
!=======================================================================
      If (irr) Then

        nsym = nsymat

!-----------------------------------------------------------------------

        If (nsym>48) Stop 'DIM problem in bzirr3d.'

        If ((lsurf) .And. (nsym>12)) Stop 'DIM problem in bzirr3d, surf mode.'

!----------------------------------------------------------------------

        Do n = 1, nsym
          isym = isymindex(n)
          msign = 1.E0_dp
          If ((krel==1) .And. (.Not. symunitary(n))) msign = -1.E0_dp
          Do i = 1, 3
            Do j = 1, 3
              u(i, j, n) = msign*rsymat(isym, i, j)
            End Do
          End Do
        End Do
!========================================================================
      Else
!========================================================================
        nsym = 1
        u(1, 1, 1) = 1.E0_dp
        u(1, 2, 1) = 0.E0_dp
        u(1, 3, 1) = 0.E0_dp
        u(2, 1, 1) = 0.E0_dp
        u(2, 2, 1) = 1.E0_dp
        u(2, 3, 1) = 0.E0_dp
        u(3, 1, 1) = 0.E0_dp
        u(3, 2, 1) = 0.E0_dp
        u(3, 3, 1) = 1.E0_dp
      End If
!==========================================================================

      Do j = 1, 3
        Do i = 1, 3
          bgmat(i, j) = ddot(3, gq(1,i), 1, gq(1,j), 1)
        End Do
      End Do

      Call rinvgj(bginv, bgmat, 3, ndim)

!-----------------------------------------------------------------------
!          rotate the 3 step vectors   GQ  --  it should be possible
!          to express the rotated vectors  B(j) in terms of the old ones
!     B(i) = SUM(j) GQ(j) * n(j,i)  with integer coefficients n(j,i)

!==========================================================================
      Do is = 1, nsym
        If (iprint>2) Write (1337, 100) is

        Do i = 1, 3
          Call dgemv('N', 3, 3, 1E0_dp, u(1,1,is), 3, gq(1,i), 1, 0E0_dp, &
            bgp(1,i), 1)

          Do j = 1, 3
            bv(j) = ddot(3, gq(1,j), 1, bgp(1,i), 1)
          End Do
          Call dgemv('N', 3, 3, 1E0_dp, bginv, 3, bv, 1, 0E0_dp, cf, 1)

          Do j = 1, ndim
            If (abs(nint(cf(j))-cf(j))>1E-8_dp) Write (1337, 120) i, j, cf(j)
            nbgp(j, i, is) = nint(cf(j))
          End Do
          If (iprint>2) Write (1337, 110) i, (bgp(j,i), j=1, 3), bv, cf, &
            (nbgp(j,i,is), j=1, 3)

        End Do

      End Do
!========================================================================

      nk = nkxyz1(1)*nkxyz1(2)*nkxyz1(3)

      Do i = 0, nkxyz1(1)
        Do j = 0, nkxyz1(2)
          Do k = 0, nkxyz1(3)
            ibk(i, j, k) = 0
          End Do
        End Do
      End Do

      ix = nkxyz1(1)
      iy = nkxyz1(2)
      iz = nkxyz1(3)
      nkp = 0
      iws = 0

      Do jx = 0, ix - 1
        iktab(1) = jx
        Do jy = 0, iy - 1
          iktab(2) = jy
          Do jz = 0, iz - 1
            iktab(3) = jz

            If (ibk(jx,jy,jz)==0) Then
              nkp = nkp + 1
              iwt = 0

!========================================================================
              Do isym = 1, nsym

!               rotate k-vector  NKP  and transform into parallelepiped
!          ROT k = SUM(i) m(i) ROT gq(i)
!                = SUM(i) m(i) SUM(j) n(i,j) gq(j)
!                = SUM(j) [SUM(i) m(i) n(i,j)] gq(j)

                Do j = 1, 3
                  is = 0
                  Do i = 1, 3
                    is = is + iktab(i)*nbgp(j, i, isym)
                  End Do
                  is = mod(is, nkxyz1(j))
                  If (is<0) is = is + nkxyz1(j)
                  ind2(j) = is
                End Do

                ja = ind2(1)
                jb = ind2(2)
                jc = ind2(3)

                If (ibk(ja,jb,jc)==0) Then
                  ibk(ja, jb, jc) = nkp
                  iwt = iwt + 1
                End If
! Mesh point in the unit cell found.
              End Do
!========================================================================

              If (nkp<=kpoibz) Then
                Do i = 1, 3
                  kp(i, nkp) = 0E0_dp
                  Do j = 1, 3
                    kp(i, nkp) = kp(i, nkp) + gq(i, j)*real(iktab(j), kind=dp)
                  End Do
                End Do
                wtkp(nkp) = real(iwt, kind=dp)/real(nk, kind=dp)
              Else
                Write (6, *) ' No. of k-points ', nkp, ' > KPOIBZ '
                Stop ' <BZIRR3D> increase parameter KPOIBZ'
              End If

            End If
            iws = iws + iwt

          End Do
        End Do
      End Do

      volbz = 0.E0_dp

      If (lsurf) Then
        v1 = abs(recbv(1,1)*recbv(2,2)-recbv(1,2)*recbv(2,1))
      Else
        v1 = ddet33(recbv)
      End If

      Do i = 1, nkp
        wtkp(i) = wtkp(i)*v1/real(nsym, kind=dp)
        volbz = volbz + wtkp(i)*real(nsym, kind=dp)
      End Do

      Deallocate (ibk)

      Return

100   Format (5X, 'rotated GQ  for IROT=', I3)
110   Format (5X, I3, 3F7.3, 2X, 3F7.3, 2X, 3F7.3, 2X, 3I3)
120   Format (5X, 2I3, 3F7.3)
    End Subroutine
