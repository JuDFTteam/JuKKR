!-------------------------------------------------------------------------------
! SUBROUTINE: lattice2d
!> @brief Generates the lattice vectors of direct and reciprocal space from
!> basic translation vectors for a 2D system
!> @note - V. Popescu May 2004: The routine has been brought to a form which is very similar to
!> LATTICE2D -- from which it has been originally derived. Dimension of arrays GN,RM
!> changed from (4,*) to (2,*), the 4th one it is used only locally (GNR/RMR) -- only GN/RM(2,*) are
!> actually needed in EWALD2D
!> @note - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
    Subroutine lattice2d(alat, bravais, recbv, ngmax, nrmax, nshlg, nshlr, &
      nsg, nsr, gn, rm, rmax, gmax, iprint, nmaxd, ishld)
! **********************************************************************
! *                                                                    *
! *  generate lattice vectors of direct and reciprocal space from      *
! *  basic translation vectors br                                      *
! *                                                                    *
! *  alat            : lattice constant                                *
! *  br(i,j)         : i=x,y,z j= 1,2,3 bravais vectors                *
! *                    *** in a.u. ****                                *
! *  rmax            : maximum radius in real space        (input)     *
! *  gmax            : maximum radius in reciprocal space  (input)     *
! *  ngmax           : Number of reciprocal lattice vectors            *
! *  gn(2,nmaxd)     : x,y,z   of reciprocal lattice vectors           *
! *  nrmax           : Number of real lattice vectors                  *
! *  rm(2,nmaxd)     : x,y,z  of real space vectors                    *
! *  nshlg           : shells in reciprocal space                      *
! *  nshlr           : shells in real space                            *
! *  nsg,nsr         : integer arrays, number of atoms in each shell   *
! *                                                                    *
! *  The routine has been brought to a form which is very similar to   *
! *  LATTICE2D -- from which it has been originally derived            *
! *  Dimension of arrays GN,RM changed from (4,*) to (2,*), the 4th    *
! *  one it is used only locally (GNR/RMR) -- only GN/RM(2,*) are      *
! *  actually needed in EWALD2D                  v.popescu May 2004    *
! *                                                                    *
! **********************************************************************
      Use mod_datatypes, Only: dp
      Implicit None
! ..
! .. Input variables
      Integer, Intent (In) :: nmaxd !< Paremeters for the Ewald summations
      Integer, Intent (In) :: ishld !< Paremeters for the Ewald summations
      Integer, Intent (In) :: iprint
      Real (Kind=dp), Intent (In) :: alat !< Lattice constant in a.u.
      Real (Kind=dp), Dimension (3, 3), Intent (In) :: recbv !< Reciprocal basis vectors
      Real (Kind=dp), Dimension (3, 3), Intent (In) :: bravais !< Bravais lattice vectors
! ..
! .. Input/Output variables
      Real (Kind=dp), Intent (Inout) :: gmax !< Ewald summation cutoff parameter for reciprocal space summation
      Real (Kind=dp), Intent (Inout) :: rmax !< Ewald summation cutoff parameter for real space summation
      Integer, Dimension (ishld), Intent (Inout) :: nsg
      Integer, Dimension (ishld), Intent (Inout) :: nsr
      Real (Kind=dp), Dimension (2, nmaxd), Intent (Inout) :: gn !< x,y,z   of reciprocal lattice vectors
      Real (Kind=dp), Dimension (2, nmaxd), Intent (Inout) :: rm !< x,y,z  of real space vectors
! .. Output variables
      Integer, Intent (Out) :: nshlr !< Shells in real space
      Integer, Intent (Out) :: nshlg !< Shells in reciprocal space
      Integer, Intent (Out) :: nrmax !< Number of real space vectors rr
      Integer, Intent (Out) :: ngmax !< Number of reciprocal space vectors
! ..
! .. Local scalars ..
      Integer :: idint
      Integer :: i, k, l, m, n, n1, ng, nr, nsh, nshl, numg, numgh, numr, &
        numrh
      Real (Kind=dp) :: dble
      Real (Kind=dp) :: rx, ry, vmin
      Real (Kind=dp) :: a, absgm, absrm, ag, ar, b, da, db, gx, gy, pi
! ..
! .. Local arrays ..
      Real (Kind=dp), Dimension (nmaxd) :: gnr
      Real (Kind=dp), Dimension (nmaxd) :: rmr
      Real (Kind=dp), Dimension (3) :: absg
      Real (Kind=dp), Dimension (3) :: absr
      Real (Kind=dp), Dimension (nmaxd) :: length
      Real (Kind=dp), Dimension (3, 3) :: bg
      Real (Kind=dp), Dimension (3, 3) :: br
      Real (Kind=dp), Dimension (4, nmaxd) :: cj
! ..
! .. Intrinsic functions ..
      Intrinsic :: abs, atan, max, mod, sqrt
! ..
! .. External subroutines ..
      External :: ioinput
!----------------------------------------------------------------------------
      pi = 4.0E0_dp*atan(1.0E0_dp)
!----------------------------------------------------------------------------
! OUTPUT
!----------------------------------------------------------------------------
      Write (1337, '(5X,2A,/)') '< LATTICE2D > : ', &
        'generating direct/reciprocal lattice vectors'
!----------------------------------------------------------------------------
! OUTPUT
!----------------------------------------------------------------------------
      rmax = rmax*alat
      gmax = gmax/alat
!----------------------------------------------------------------------------
! OUTPUT
!----------------------------------------------------------------------------
      Write (1337, Fmt=100) rmax, gmax
!----------------------------------------------------------------------------
! OUTPUT
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
! Basic trans. vectors and basis vectors
!----------------------------------------------------------------------------
      Do i = 1, 3
        br(1, i) = bravais(1, i)*alat
        br(2, i) = bravais(2, i)*alat
      End Do
!----------------------------------------------------------------------------
! Generate primitive vectors BG of reciprocal space
!----------------------------------------------------------------------------
      Do i = 1, 3
        bg(1, i) = recbv(1, i)*2E0_dp*pi/alat
        bg(2, i) = recbv(2, i)*2E0_dp*pi/alat
      End Do
!----------------------------------------------------------------------------
! Estimate no. of lattice vectors
!----------------------------------------------------------------------------
      Do i = 1, 3
        absr(i) = sqrt(br(1,i)**2+br(2,i)**2)
        absg(i) = sqrt(bg(1,i)**2+bg(2,i)**2)
      End Do
!
      absrm = max(absr(1), absr(2))
      absgm = max(absg(1), absg(2))
      absrm = 2.0E0_dp*pi/absrm
      absgm = 2.0E0_dp*pi/absgm
      numr = 2*(idint(rmax/absgm)+1) + 1
      numg = 2*(idint(gmax/absrm)+1) + 1
      numrh = numr/2 + 1
      numgh = numg/2 + 1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 generate lattice vectors of real space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      Write (1337, *) 'Real space...'
      nr = 0
!----------------------------------------------------------------------------
      Do l = 1, numr
        a = dble(l-numrh)
        Do m = 1, numr
          b = dble(m-numrh)
!----------------------------------------------------------------------
          rx = a*br(1, 1) + b*br(1, 2)
          ry = a*br(2, 1) + b*br(2, 2)
          ar = sqrt(rx*rx+ry*ry)
!----------------------------------------------------------------------
          If (ar<=rmax) Then
            nr = nr + 1
            If (nr>nmaxd) Then
              Write (6, *) &
                'lattice2d: ERROR: Dimension NMAXD in the inputcard too small' &
                , nr, nmaxd
              Stop 'lattice2d'
            End If
            cj(1, nr) = rx
            cj(2, nr) = ry
            cj(3, nr) = 0E0_dp
            cj(4, nr) = ar
          End If
        End Do
      End Do
!----------------------------------------------------------------------------
      nrmax = nr
!----------------------------------------------------------------------------
! Sort vectors in order of increasing absolute value
!----------------------------------------------------------------------------
      Write (1337, Fmt='(A11,I8,A11)') '...sorting ', nrmax, ' vectors...'

      da = 1.E-06_dp
      nsh = 0
      nshl = -1
!----------------------------------------------------------------------------
      Do k = 1, nr
        vmin = rmax + 1.0E0_dp
        Do n = 1, nr
          If (cj(4,n)-vmin<0E0_dp) Then
            vmin = cj(4, n)
            n1 = n
          End If
        End Do
!
        nshl = nshl + 1
        rm(1, k) = cj(1, n1)
        rm(2, k) = cj(2, n1)
        rmr(k) = cj(4, n1)
        db = vmin
!-------------------------------------------------------------------------
        If (db>da+1.E-06_dp) Then
          nsh = nsh + 1
          If (nsh>ishld) Then
            Write (6, *) ' ERROR: Dimension ISHLD in the inputcard too small', &
              nsh, ishld
            Stop 'lattice2d'
          End If
!
          nsr(nsh) = nshl
          nshl = 0
          da = db
        End If
!-------------------------------------------------------------------------
        cj(4, n1) = rmax + 1.0E0_dp
      End Do
!----------------------------------------------------------------------------
      nsh = nsh + 1
      nshl = nshl + 1
      If (nsh>ishld) Then
        Write (6, *) ' ERROR: Dimension ISHLD in the inputcard too small', &
          nsh, ishld
        Stop 'lattice2d'
      End If
!
      nsr(nsh) = nshl
      nshlr = nsh
      If (nshlr<=1) Stop 'lattice2d: ERROR: cut-off radius RMAX too small '
      Write (1337, *) '...done.'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 generate lattice vectors of real space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      Write (1337, *) 'Reciprocal space...'
      ng = 0
!----------------------------------------------------------------------------
      Do l = 1, numg
        a = dble(l-numgh)
        Do m = 1, numg
          b = dble(m-numgh)
!----------------------------------------------------------------------
          gx = a*bg(1, 1) + b*bg(1, 2)
          gy = a*bg(2, 1) + b*bg(2, 2)
          ag = sqrt(gx*gx+gy*gy)
!----------------------------------------------------------------------
          If (ag<=gmax) Then
            ng = ng + 1
            If (ng>nmaxd) Then
              Write (6, *) &
                ' ERROR: Dimension NMAXD in the inputcard too small', ng, &
                nmaxd
              Stop 'lattice2d'
            End If
            cj(1, ng) = gx
            cj(2, ng) = gy
            cj(3, ng) = 0E0_dp
            cj(4, ng) = ag
          End If
        End Do
      End Do
!----------------------------------------------------------------------------
      ngmax = ng
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
! Sort vectors in order of increasing abs. value
!----------------------------------------------------------------------------
      Write (1337, Fmt='(A11,I8,A11)') '...sorting ', ngmax, ' vectors...'
      Do n = 1, ng
        length(n) = cj(4, n)
      End Do
      da = 1.E-06_dp
      nsh = 0
      nshl = -1
!----------------------------------------------------------------------------
      Do k = 1, ng
        vmin = gmax + 1.0E0_dp
        Do n = 1, ng
          If (length(n)<vmin) Then ! ( CJ(4,N).LT.VMIN ) THEN
            vmin = length(n) ! CJ(4,N)
            n1 = n
          End If
        End Do
!
        nshl = nshl + 1
        gn(1, k) = cj(1, n1)
        gn(2, k) = cj(2, n1)
        gnr(k) = length(n1) ! CJ(4,N1)
        db = vmin
!-------------------------------------------------------------------------
        If (db>da+1.E-07_dp) Then
          nsh = nsh + 1 ! Number of shells of different length
          If (nsh>ishld) Then
            Write (6, *) ' ERROR: Dimension ISHLD in the inputcard too small', &
              nsh, ishld
            Stop 'lattice2d'
          End If
!
          nsg(nsh) = nshl ! Number of vectors in shell
          nshl = 0
          da = db
        End If
!-------------------------------------------------------------------------
        length(n1) = gmax + 1.0E0_dp !  CJ(4,N1) = GMAX + 1.0D0
      End Do
!----------------------------------------------------------------------------
      nsh = nsh + 1
      nshl = nshl + 1
      If (nsh>ishld) Then
        Write (6, *) ' ERROR: Dimension ISHLD in the inputcard too small', &
          nsh, ishld
        Stop 'lattice2d'
      End If
!
      nsg(nsh) = nshl
      nshlg = nsh
      If (nshlg<=1) Stop 'lattice2dERROR: cut-off radius GMAX too small '

      Write (1337, *) '...done.'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! OUTPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Write (1337, Fmt=110)
      Write (1337, Fmt=120) 'Direct  lattice', nrmax, nshlr, rmr(nrmax)
      Write (1337, Fmt=120) 'Recipr. lattice', ngmax, nshlg, gnr(ngmax)
      Write (1337, Fmt=130)
!
      If (iprint<3) Return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! OUTPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------------
      k = 0
      Write (1337, Fmt=140) 'real-space'
      Do l = 1, nshlr
        Write (1337, 150) l, nsr(l), rmr(k+1), (rm(m,k+1), m=1, 2)
        Do n = 2, nsr(l)
          Write (1337, Fmt=160)(rm(m,k+n), m=1, 2)
        End Do
        If (l/=nshlr) Write (1337, 170)
        k = k + nsr(l)
      End Do
      Write (1337, 180)
      k = 0
      Write (1337, Fmt=140) 'reciprocal'
      Do l = 1, nshlg
        Write (1337, 150) l, nsg(l), gnr(k+1), (gn(m,k+1), m=1, 2)
        Do n = 2, nsg(l)
          Write (1337, Fmt=160)(gn(m,k+n), m=1, 2)
        End Do
        If (l/=nshlg) Write (1337, 170)
        k = k + nsg(l)
      End Do
      Write (1337, 180)
!----------------------------------------------------------------------------
!
100   Format (10X, 'R max =', F10.5, ' (a.u.)', /, 10X, 'G max =', F10.5, &
        ' (1/a.u.)', /)
110   Format (10X, '               vectors  shells  max. R ', /, 10X, &
        '               ------------------------------')
120   Format (10X, A, I7, 2X, I6, 2X, F9.5)
130   Format (10X, '               ------------------------------', /)
140   Format (10X, 45('+'), /, 13X, 'generated ', A, ' lattice vectors', /, &
        10X, 45('+'), /, 10X, 'shell Nvec    radius          x         y', /, &
        10X, 45('-'))
150   Format (10X, I5, I5, F12.6, 2X, 2F10.5)
160   Format (34X, 2F10.5)
170   Format (13X, 42('-'))
180   Format (10X, 45('+'), /)
    End Subroutine
