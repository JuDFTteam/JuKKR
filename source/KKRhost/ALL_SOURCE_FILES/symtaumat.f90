    Subroutine symtaumat(rotname, rotmat, drot, nsym, isymindex, symunitary, &
      nqmax, nkmmax, nq, nl, krel, iprint, nsymaxd)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *  Find the symmetry matrices DROT that act on t, tau, ....        *
!   *  KREL=0: for real spherical harmonics                            *
!   *  KREL=1: for relativistic represntation                          *
!   *                                                                  *
!   *  The NSYM allowed symmetry operations are indicated by ISYMINDEX *
!   *  in the table  ROTMAT. For KREL=1, SYMUNITARY=T/F indicates a    *
!   *  unitary/antiunitary symmetry operation.                         *
!   *                                                                  *
!   *  The routine determines first the Euler angles correponding      *
!   *  to a symmetry operation. Reflections are decomposed into        *
!   *  inversion + rotation for this reason.                           *
!   *                                                                  *
!   ********************************************************************

      Implicit None

! PARAMETER definitions
      Complex (Kind=dp) :: ci, c1, c0
      Parameter (ci=(0.0E0_dp,1.0E0_dp), c1=(1.0E0_dp,0.0E0_dp), &
        c0=(0.0E0_dp,0.0E0_dp))

! Dummy arguments
      Integer :: iprint, krel, nkmmax, nl, nq, nqmax, nsym, nsymaxd
      Complex (Kind=dp) :: drot(nkmmax, nkmmax, 48)
      Integer :: isymindex(nsymaxd)
      Real (Kind=dp) :: rotmat(64, 3, 3)
      Character (Len=10) :: rotname(64)
      Logical :: symunitary(48)

! Local variables
      Real (Kind=dp) :: a, b, co1, co2, co3, det, fact(0:100), pi, si1, si2, &
        si3, sk, symeulang(3, 48), tet1, tet2, tet3
      Logical :: checkrmat
      Real (Kind=dp) :: dble
      Real (Kind=dp) :: ddet33
      Complex (Kind=dp) :: dinv(nkmmax, nkmmax), dtim(nkmmax, nkmmax), &
        rc(nkmmax, nkmmax), w1(nkmmax, nkmmax), w2(nkmmax, nkmmax)
      Logical :: equal
      Integer :: i, i1, i2, ind0q(nqmax), invflag(48), iq, irel, ireleff, &
        isym, itop, j, k, l, loop, m, n, nk, nkeff, nkm, nlm, nok, ns, rj, rmj
      Integer :: nint
      Real (Kind=dp) :: rmat(3, 3)
      Real (Kind=dp) :: w

      equal(a, b) = (dabs(a-b)<1E-7_dp)

      Write (1337, 110)

      pi = 4E0_dp*atan(1E0_dp)

      irel = krel*3
      nk = (1-krel)*nl + krel*(2*nl-1)
      nkm = (1+krel)*nl**2

!-----------------------------------------------------------------------
      fact(0) = 1.0E0_dp
      Do i = 1, 100
        fact(i) = fact(i-1)*dble(i)
      End Do
!-----------------------------------------------------------------------

      ind0q(1) = 0
      Do iq = 2, nq
        ind0q(iq) = ind0q(iq-1) + nkm
      End Do

! ----------------------------------------------------------------------
!    RC  transforms from  REAL to  COMPLEX (L,M,S) - representation
!                 |LC> = sum[LR] |LR> * RC(LR,LC)
! ----------------------------------------------------------------------
      If (krel==0) Then
        nlm = nkm

        Call cinit(nkmmax*nkmmax, rc)

        w = 1.0E0_dp/sqrt(2.0E0_dp)

        Do l = 0, (nl-1)
          Do m = -l, l
            i = l*(l+1) + m + 1
            j = l*(l+1) - m + 1

            If (m<0) Then
              rc(i, i) = -ci*w
              rc(j, i) = w
            End If
            If (m==0) Then
              rc(i, i) = c1
            End If
            If (m>0) Then
              rc(i, i) = w*(-1.0E0_dp)**m
              rc(j, i) = ci*w*(-1.0E0_dp)**m
            End If
          End Do
        End Do
      End If

!=======================================================================
!     The routine determines first the Euler angles correponding
!     to a symmetry operation. Reflections are decomposed into
!     inversion + rotation for this reason.
!=======================================================================

      Do isym = 1, nsym

        Do i1 = 1, 3
          Do i2 = 1, 3
            rmat(i1, i2) = rotmat(isymindex(isym), i1, i2)
          End Do
        End Do

        det = ddet33(rmat)

        invflag(isym) = 0
        If (det<0E0_dp) Then
          Call dscal(9, -1.0E0_dp, rmat, 1)
          invflag(isym) = 1
        End If

!----------------------------------------------------------------------
        co2 = rmat(3, 3)
        tet2 = acos(co2)
        loop = 0
100     Continue
        If (loop==1) tet2 = -tet2
        si2 = sin(tet2)

        If (equal(co2,1.0E0_dp)) Then
          tet1 = acos(rmat(1,1))
          If (.Not. equal(rmat(1,2),sin(tet1))) Then
            tet1 = -tet1
            If (.Not. equal(rmat(1,2),sin(tet1))) Write (1337, *) &
              '>>>>>>>>>>>>>>> STRANGE 1'
          End If
          tet2 = 0E0_dp
          tet3 = 0E0_dp
        Else If (equal(co2,-1E0_dp)) Then
          tet1 = acos(-rmat(1,1))
          If (.Not. equal(rmat(1,2),-sin(tet1))) Then
            tet1 = -tet1
            If (.Not. equal(rmat(1,2),-sin(tet1))) Write (1337, *) &
              '>>>>>>>>>>>>>>> STRANGE 2'
          End If
          tet2 = pi
          tet3 = 0E0_dp
        Else
          tet1 = acos(rmat(3,1)/si2)
          If (.Not. equal(rmat(3,2),si2*sin(tet1))) Then
            tet1 = -tet1
            If (.Not. equal(rmat(3,2),si2*sin(tet1))) Write (1337, *) &
              '>>>>>>>>>>>>>>> STRANGE 3'
          End If

          tet3 = acos(-rmat(1,3)/si2)
          If (.Not. equal(rmat(2,3),si2*sin(tet3))) Then
            tet3 = -tet3
            If (.Not. equal(rmat(2,3),si2*sin(tet3))) Write (1337, *) &
              '>>>>>>>>>>>>>>> STRANGE 4'
          End If

        End If

        co1 = cos(tet1)
        si1 = sin(tet1)
        co2 = cos(tet2)
        si2 = sin(tet2)
        co3 = cos(tet3)
        si3 = sin(tet3)

        nok = 0
        Do i1 = 1, 3
          Do i2 = 1, 3
            If (checkrmat(rmat,co1,si1,co2,si2,co3,si3,i1,i2)) Then
              nok = nok + 1
            Else If (loop<1) Then
              loop = loop + 1
              Go To 100
            End If
          End Do
        End Do

        symeulang(1, isym) = tet1*(180E0_dp/pi)
        symeulang(2, isym) = tet2*(180E0_dp/pi)
        symeulang(3, isym) = tet3*(180E0_dp/pi)

        If (nok/=9) Write (1337, 130) nok
        Write (1337, 120) isym, rotname(isymindex(isym)), invflag(isym), &
          (symeulang(i,isym), i=1, 3), symunitary(isym)

      End Do
      Write (1337, '(8X,57("-"),/)')

!-----------------------------------------------------------------------
!                    initialize all rotation matrices
!-----------------------------------------------------------------------

      Call cinit(nkmmax*nkmmax*nsym, drot)

!-----------------------------------------------------------------------
!                       create rotation matrices
!-----------------------------------------------------------------------

      If (irel<=2) Then
        ireleff = 0
        nkeff = nl
      Else
        ireleff = 3
        nkeff = nk
      End If

      Do isym = 1, nsym

        Call calcrotmat(nkeff, ireleff, symeulang(1,isym), symeulang(2,isym), &
          symeulang(3,isym), drot(1,1,isym), fact, nkmmax)

      End Do
!-----------------------------------------------------------------------
!                     create matrix for inversion
!-----------------------------------------------------------------------
      Call cinit(nkmmax*nkmmax, dinv)

      i = 0
      If (irel>2) Then
        ns = 2
      Else
        ns = 1
      End If
      Do l = 0, (nl-1)
        Do m = 1, ns*(2*l+1)
          i = i + 1
          dinv(i, i) = (-1.0E0_dp)**l
        End Do
      End Do
      itop = i

!-----------------------------------------------------------------------
!                         include inversion
!-----------------------------------------------------------------------
      Do isym = 1, nsym
        If (invflag(isym)/=0) Then

          Call zgemm('N', 'N', nkm, nkm, nkm, c1, drot(1,1,isym), nkmmax, &
            dinv, nkmmax, c0, w2, nkmmax)

          Do j = 1, nkm
            Call zcopy(nkm, w2(1,j), 1, drot(1,j,isym), 1)
          End Do
        End If
      End Do

!-----------------------------------------------------------------------
!            add second spin-diagonal block for  IREL=2
!            spin off-diagonal blocks have been initialized before
!-----------------------------------------------------------------------
      If (irel==2) Then
        nlm = nkm/2
        If (itop/=nlm) Call errortrap('SYMTAUMAT', 11, 1)
        Do isym = 1, nsym

          Do j = 1, nlm
            Call zcopy(nlm, drot(1,j,isym), 1, drot(nlm+1,nlm+j,isym), 1)
          End Do
        End Do
      End If
!-----------------------------------------------------------------------
!            transform to real spherical representation for  KREL=0
!-----------------------------------------------------------------------
      n = nkm
      m = nkmmax
      If (krel==0) Then
        Do isym = 1, nsym
          Call zgemm('N', 'N', n, n, n, c1, rc, m, drot(1,1,isym), m, c0, w1, &
            m)
          Call zgemm('N', 'C', n, n, n, c1, w1, m, rc, m, c0, drot(1,1,isym), &
            m)
        End Do
      End If
!-----------------------------------------------------------------------
!                     create matrix for time reversal
!-----------------------------------------------------------------------
      If (irel>1) Then

        Call cinit(nkmmax*nkmmax, dtim)

        i = 0
        Do k = 1, nk
          l = k/2
          If (l*2==k) Then
            sk = -1E0_dp
          Else
            sk = +1E0_dp
          End If
          rj = nint(l+sk*0.5E0_dp)
          Do rmj = -rj, +rj, 1
            i1 = nint(2*l*(rj+0.5E0_dp)+rj+rmj+1)
            i2 = nint(2*l*(rj+0.5E0_dp)+rj-rmj+1)
            dtim(i1, i2) = sk*(-1)**nint(rmj+0.5E0_dp)
          End Do
        End Do
        If (iprint>0) Then
          Call cmatstr('Inversion     MATRIX', 20, dinv, nkm, nkmmax, 3, 3, 0, &
            1E-8_dp, 6)
          Call cmatstr('Time reversal MATRIX', 20, dtim, nkm, nkmmax, 3, 3, 0, &
            1E-8_dp, 6)
        End If

      End If
!=======================================================================
!            set up of transformation matrices completed
!=======================================================================

!=======================================================================
!   include time reversal operation for anti-unitary transformations
!=======================================================================
      Do isym = 1, nsym
        If (.Not. symunitary(isym)) Then
          If (irel==2) Call errortrap('SYMTAUMAT', 14, 1)

          Call zgemm('N', 'N', nkm, nkm, nkm, c1, drot(1,1,isym), nkmmax, &
            dtim, nkmmax, c0, w2, nkmmax)
          Do j = 1, nkm
            Call zcopy(nkm, w2(1,j), 1, drot(1,j,isym), 1)
          End Do
        End If
      End Do

!-----------------------------------------------------------------------
! for testing

!ccc      write (6,*) ' NUMBER OF SYMMETRIES : ', NSYM
!ccc
!ccc      do isym = 1,nsym
!ccc         write(6,*) ' ISYM = ',isym
!ccc         call cmatstr('DROT',4,drot(1,1,isym),nkm,nkmmax,krel*3,krel*3,
!ccc     &        0,1d-12,6)
!ccc         write(6,*)
!ccc      end do

!-----------------------------------------------------------------------

      If (iprint==0) Return

!=======================================================================
!       find the structure of the site-diagonal TAU - matrices  TAUQ
!=======================================================================

      Call taustruct(drot, nsym, symunitary, nkm, nq, nqmax, nkmmax, iprint, &
        irel)

      Return
110   Format (5X, '<SYMTAUMAT> : rotation matrices acting on t/G/tau', /, /, &
        8X, 57('-'), /, 8X, &
        'ISYM            INV          Euler angles      Unitarity', /, 8X, &
        57('-'))
120   Format (8X, I2, 3X, A, I3, 3F10.5, 3X, L1)
130   Format (50('>'), ' trouble in <SYMTAUMAT>', I3, F10.5)
    End Subroutine
