    Subroutine mvecglobal(it, iq, natyp, qmphi, qmtet, mvevi, mvevil, mvevief, &
      natypd, lmaxd, nmvecmax)
!   ********************************************************************
!   *                                                                  *
!   *  this routine has been build up from the last part of the        *
!   *  original Munich CALCMVEC routine.                               *
!   *  on exit, MVEVI,MVEVIL,MVEVIEF are in the CARTESIAN GLOBAL       *
!   *                                           frame of reference     *
!   *                                                                  *
!   ********************************************************************
      Use mod_datatypes
      Implicit None

!Parameter definitions
      Integer :: lmaxdloc
      Parameter (lmaxdloc=8)
      Complex (Kind=dp) :: ci, czero
      Parameter (ci=(0.0E0_dp,1.0E0_dp), czero=(0.0E0_dp,0.0E0_dp))

!Scalar Arguments
      Integer :: it, iq, natyp, natypd, lmaxd, nmvecmax
      Real (Kind=dp) :: qmphi, qmtet

!Array Arguments
      Complex (Kind=dp) :: mvevi(natypd, 3, nmvecmax), &
        mvevil(0:lmaxd, natypd, 3, nmvecmax)
      Complex (Kind=dp) :: mvevief(natypd, 3, nmvecmax)

!Local Scalars
      Integer :: icall, i, j, k, l, imv, nmvec
      Complex (Kind=dp) :: cs
      Complex (Kind=dp) :: amin, apls
      Real (Kind=dp) :: pi, mv, mvx, mvxy, mvy, mvz, wsq2

!Local Arrays
      Complex (Kind=dp) :: usc(3, 3), drot4(4, 4), w3x3(3, 3)
      Complex (Kind=dp) :: mvg(3, nmvecmax), mvgef(3, nmvecmax)
      Complex (Kind=dp) :: mvgl(0:lmaxd, 3, nmvecmax)
      Real (Kind=dp) :: mrot(3, 3), fact(0:100)
      Real (Kind=dp) :: mvglo(3, nmvecmax), mvglol(0:lmaxd, 3, nmvecmax)
      Real (Kind=dp) :: mvphi(nmvecmax), mvtet(nmvecmax)
      Character (Len=1) :: txtl(0:lmaxdloc)

!Intrinsic Functions
      Intrinsic :: abs, acos, atan, real, aimag, conjg

!External subroutines
      External :: calcrotmat

!Data Statements
      Data icall/0/

!Save Statements
      Save :: icall, nmvec, usc, fact, txtl, pi

      icall = icall + 1
!=======================================================================
      If (icall==1) Then

        If (lmaxd>lmaxdloc) Then
          Write (6, *)
          Write (6, *) ' Please increase parameter LMAXDLOC to ', lmaxd
          Write (6, *) ' in the < MVECGLOBAL > routine.'
          Stop ' < TBKKR2 > '
        End If

        txtl(0) = 's'
        txtl(1) = 'p'
        txtl(2) = 'd'
        If (lmaxd>=3) Then
          Do l = 3, lmaxd
            txtl(l) = char(ichar('f')+l-3)
          End Do
        End If

        Write (1337, '(78("#"))')
        Write (1337, 100)
        Write (1337, '(78("#"))')
        Write (1337, *)
        Write (1337, 110)

        nmvec = 2
        pi = 4.E0_dp*atan(1.E0_dp)

        fact(0) = 1.0E0_dp
        Do i = 1, 100
          fact(i) = fact(i-1)*real(i, kind=dp)
        End Do
!-----------------------------------------------------------------------
!  create transformation matrix   U  cartesian/sperical ccordinates
!-----------------------------------------------------------------------
!  RC,RCP  vectors in cartesian coordinates
!  RS,RSP  vectors in spherical coordinates
!         RS  = USC * R!..                          (4.40)
!         RSP = MS  * RS                                 (4.37)
!     MS(i,j) = D(j,i)                                   (4.42)
!     D  rotation matrix for complex spherical harmonics

! ordering of: m=-1,0,+1 >>> row 1 and 3 interchanged compared to (4.44)

        wsq2 = 1.0E0_dp/sqrt(2.0E0_dp)

        usc(1, 1) = wsq2
        usc(1, 2) = -ci*wsq2
        usc(1, 3) = 0.0E0_dp
        usc(2, 1) = 0.0E0_dp
        usc(2, 2) = 0.0E0_dp
        usc(2, 3) = 1.0E0_dp
        usc(3, 1) = -wsq2
        usc(3, 2) = -ci*wsq2
        usc(3, 3) = 0.0E0_dp
!-----------------------------------------------------------------------
      End If
!=======================================================================

!-----------------------------------------------------------------------
!   create the rotation matrices  DROT4 for complex spherical harmonics
!-----------------------------------------------------------------------

      Call calcrotmat(2, 1, qmphi, qmtet, 0.0E0_dp, drot4, fact, 4)

!-----------------------------------------------------------------------
! create the rotation matrix  MROT for vectors in cartesian coordinates
! NOTE:  U^+ D^T U gives the inverse of the real matrix  M
!        for that reason  the transposed matrix is stored as  MROT(J,I)
!-----------------------------------------------------------------------

      Do i = 1, 3
        Do j = 1, 3
          cs = 0.0E0_dp
          Do k = 1, 3
            cs = cs + drot4(k+1, i+1)*usc(k, j)
          End Do
          w3x3(i, j) = cs
        End Do
      End Do

      Do i = 1, 3
        Do j = 1, 3
          cs = 0.0E0_dp
          Do k = 1, 3
            cs = cs + conjg(usc(k,i))*w3x3(k, j)
          End Do
          If (aimag(cs)>1E-8_dp) Write (*, *) ' MROT', i, j, cs, &
            ' ???????????'
!     see above >> MROT(I,J) = DREAL(CS)
          mrot(j, i) = real(cs)
        End Do
      End Do
!-----------------------------------------------------------------------

! **********************************************************************
      Do imv = 1, nmvec
!-----------------------------------------------------------------------
!     transform from (+,-,z) to cartesian coordinates  (x,y,z)
!     note the convention
!-----------------------------------------------------------------------
        apls = mvevi(it, 1, imv)
        amin = mvevi(it, 2, imv)
        mvevi(it, 1, imv) = (amin+apls)*0.5E0_dp
        mvevi(it, 2, imv) = (amin-apls)*0.5E0_dp*ci

        apls = mvevief(it, 1, imv)
        amin = mvevief(it, 2, imv)
        mvevief(it, 1, imv) = (amin+apls)*0.5E0_dp
        mvevief(it, 2, imv) = (amin-apls)*0.5E0_dp*ci

        Do l = 0, lmaxd
          apls = mvevil(l, it, 1, imv)
          amin = mvevil(l, it, 2, imv)
          mvevil(l, it, 1, imv) = (amin+apls)*0.5E0_dp
          mvevil(l, it, 2, imv) = (amin-apls)*0.5E0_dp*ci
        End Do
!-----------------------------------------------------------------------
!     transform from LOCAL cartesian coordinates (x,y,z)
!               to  GLOBAL cartesian coordinates
!-----------------------------------------------------------------------
        Do i = 1, 3
          mvg(i, imv) = czero
          mvgef(i, imv) = czero
          Do j = 1, 3
            mvg(i, imv) = mvg(i, imv) + mrot(i, j)*mvevi(it, j, imv)
            mvgef(i, imv) = mvgef(i, imv) + mrot(i, j)*mvevief(it, j, imv)
          End Do
          mvglo(i, imv) = aimag(mvg(i,imv))

          Do l = 0, lmaxd
            mvgl(l, i, imv) = czero
            Do j = 1, 3
              mvgl(l, i, imv) = mvgl(l, i, imv) + mrot(i, j)*mvevil(l, it, j, &
                imv)
            End Do
            mvglol(l, i, imv) = aimag(mvgl(l,i,imv))
          End Do

        End Do
! ......................................................................
        Do i = 1, 3
          mvevi(it, i, imv) = mvg(i, imv)
          mvevief(it, i, imv) = mvgef(i, imv)

          Do l = 0, lmaxd
            mvevil(l, it, i, imv) = mvgl(l, i, imv)
          End Do
        End Do
!-----------------------------------------------------------------------
!        calculate the angles
!-----------------------------------------------------------------------
        mvx = mvglo(1, imv)
        mvy = mvglo(2, imv)
        mvz = mvglo(3, imv)

        mv = sqrt(mvx**2+mvy**2+mvz**2)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        If (mv<1E-8_dp) Then
          mvphi(imv) = 0E0_dp
          mvtet(imv) = 0E0_dp
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        Else
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          mvxy = sqrt(mvx**2+mvy**2)
! ======================================================================
          If (abs(mvxy)<1E-8_dp) Then
            mvphi(imv) = 0E0_dp
! ======================================================================
          Else
! ======================================================================
            If (mvy>=0E0_dp) Then
              mvphi(imv) = acos(mvx/mvxy)
            Else If (mvx<0E0_dp) Then
              mvphi(imv) = pi + acos(-mvx/mvxy)
            Else
              mvphi(imv) = 2*pi - acos(mvx/mvxy)
            End If
            mvphi(imv) = mvphi(imv)*180E0_dp/pi
            If (abs(mvphi(imv)-360.0E0_dp)<1E-8_dp) mvphi(imv) = 0E0_dp
          End If
! ======================================================================
          If (mvphi(imv)>=345.E0_dp) mvphi(imv) = 360.E0_dp - mvphi(imv)
          mvtet(imv) = acos(mvz/mv)*180E0_dp/pi
        End If
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      End Do
! **********************************************************************
! output vector components, in and out angles
! ----------------------------------------------------------------------
      l = 0
      Write (1337, 120) it, iq, txtl(l), ((mvglol(l,i,imv),i=1,3), imv=1, 2)
      Write (1337, 130)(txtl(l), ((mvglol(l,i,imv),i=1, &
        3),imv=1,2), l=1, lmaxd)
      Write (1337, 140)((mvglo(i,imv),i=1,3), imv=1, 2)
      Write (1337, 150) qmphi, qmtet, (mvphi(imv), mvtet(imv), imv=1, 2)
! ----------------------------------------------------------------------
      If (it<natyp) Then
        Write (1337, '(3X,75("="))')
      Else
        Write (1337, *)
        Write (1337, '(78("#"))')
      End If
! ----------------------------------------------------------------------

100   Format (15X, 'vectorial magnetic properties given with respect', /, 15X, &
        '   to the GLOBAL (crystal) frame of reference')
110   Format (29X, 'm_spin', 27X, 'm_orb', /, 3X, 'ATOM/SITE     ', &
        '    x         y         z', 8X, '    x         y         z', /, 3X, &
        75('='))
120   Format (3X, I3, '/', I3, 2X, A1, ' =', 3F10.5, 3X, 3F10.5)
130   Format (12X, A1, ' =', 3F10.5, 3X, 3F10.5)
140   Format (12X, 66('-'), /, 12X, 'sum', 3F10.5, 3X, 3F10.5, /)
150   Format (3X, 'angles (IN)   TET =', F9.4, ' PHI =', F9.4, /, 3X, &
        'angles (calc) TET =', F9.4, ' PHI =', F9.4, '   TET =', F9.4, &
        ' PHI =', F9.4)
    End Subroutine
