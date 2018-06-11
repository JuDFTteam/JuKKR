! 02.08.95 *************************************************************
    Subroutine rrgen(bv1, lsurf, rr, nr, nrd)
      Use mod_datatypes, Only: dp
! **********************************************************************
! *                                                                    *
! * generates a number of real space vectors to construct the          *
! * clusters representing the local surrounding of the atoms in        *
! * routine CLSGEN99                                                   *
! *                                                                    *
! **********************************************************************
      Implicit None
!..
!.. Scalar arguments ..
      Logical :: lsurf
      Integer :: nr, nrd
!    ..
!    .. Array arguments ..
      Real (Kind=dp) :: bv1(3, 3), rr(3, 0:nrd)
!    ..
!    .. Local scalars ..
      Real (Kind=dp) :: epsshl, r, r1, r2, r3, rmax, rr2, rs
      Integer :: i, j, k, n1, n2, n3, pos, iprint
      Integer :: nint
      Real (Kind=dp) :: dble
!..
!.. Local arrays
      Real (Kind=dp) :: rabs(nrd), rr1(3, nrd), v(3), vx(3), vy(3), vz(3), &
        vx0(3), vy0(3), vz0(3)
      Integer :: ind(nrd)
!..
!.. Intrinsic Functions ..
      Intrinsic :: abs, min, sqrt, nint
!..
!.. External Subroutines ..
      External :: dsort, scalpr, vadd, veq
!..
!.. Data Statements ..
      Data epsshl/1.0E-5_dp/
!     ..................................................................
      Write (1337, '(5X,A,/)') &
        '< RRGEN > : generation of real space mesh RR(NR)'

      iprint = 0

      Call scalpr(bv1(1,1), bv1(1,1), r1)
      Call scalpr(bv1(1,2), bv1(1,2), r2)
      Call scalpr(bv1(1,3), bv1(1,3), r3)
      rmax = 5.E0_dp

      r1 = sqrt(r1)
      r2 = sqrt(r2)
      r3 = sqrt(r3)
      r = 1.5E0_dp*rmax + sqrt(r1*r1+r2*r2+r3*r3) + epsshl
      rs = r*r
      n1 = nint(r/r1)
      n2 = nint(r/r2)
      If (.Not. lsurf) n3 = nint(r/r3)

      n1 = min(12, n1)
      n2 = min(12, n2)
      If (.Not. lsurf) n3 = min(12, n3)

      n1 = max(2, n1)
      n2 = max(2, n2)
      If (.Not. lsurf) n3 = max(2, n3)

      If (lsurf) n3 = 0

      Write (1337, 100) r
      Write (1337, 110) rs
      If (lsurf) Then
        Write (1337, 120) n1, n2
      Else
        Write (1337, 130) n1, n2, n3
      End If

      nr = 0
      rr(1, 0) = 0.0E0_dp
      rr(2, 0) = 0.0E0_dp
      rr(3, 0) = 0.0E0_dp

      Call vmul(bv1(1,1), dble(-n1-1), vx0(1))
      Call vmul(bv1(1,2), dble(-n2-1), vy0(1))
      Call vmul(bv1(1,3), dble(-n3-1), vz0(1))
      Call veq(vx0, vx)
! **********************************************************************
      Do i = -n1, n1
        Call vadd(vx, bv1(1,1), vx)
        Call veq(vy0, vy)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        Do j = -n2, n2
          Call vadd(vy, bv1(1,2), vy)
          Call veq(vz0, vz)
! ----------------------------------------------------------------------
          Do k = -n3, n3
            Call vadd(vz, bv1(1,3), vz)
            Call vadd(vx, vy, v)
            Call vadd(v, vz, v)
            Call scalpr(v, v, rr2)

            If (((rr2<=rs) .Or. (abs(i)+abs(j)+abs(k)<= &
              6)) .And. (rr2>epsshl)) Then
              nr = nr + 1

              If (nr>nrd) Then
                Write (6, *) 'Dimension ERROR. Please, change the ', &
                  'parameter NRD in inc.p to ', nr, nrd
                Stop
              End If

              rr1(1, nr) = v(1)
              rr1(2, nr) = v(2)
              rr1(3, nr) = v(3)
              rabs(nr) = sqrt(rr2)
            End If
          End Do
! ----------------------------------------------------------------------
        End Do
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      End Do
! **********************************************************************

      Write (1337, 140) nr + 1

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      If (iprint>0) Then
        Write (1337, 150)
        Write (1337, 170) 0, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp
      End If
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

      Call dsort(rabs, ind, nr, pos)
      Do i = 1, nr
        pos = ind(i)
        rr(1, i) = rr1(1, pos)
        rr(2, i) = rr1(2, pos)
        rr(3, i) = rr1(3, pos)
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
        If (iprint>0) Write (1337, 170) i, rr(1, i), rr(2, i), rr(3, i), &
          rabs(pos)
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      End Do

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      If (iprint>0) Write (1337, 160)
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

100   Format (10X, 'Radius R        : ', F15.6, ' (ALAT    units)')
110   Format (10X, '       R**2     : ', F15.6, ' (ALAT**2 units)')
120   Format (10X, 'mesh divisions  : ', 5X, 2I5)
130   Format (10X, 'mesh divisions  : ', 3I5)
140   Format (10X, 'vectors created : ', I15)
150   Format (/, 10X, 60('+'), /, 18X, &
        'generated real-space mesh-points (ALAT units)', /, 10X, 60('+'), /, &
        13X, 'index      x           y           z          distance  ', /, &
        10X, 60('-'))
160   Format (10X, 60('+'))
170   Format (10X, I6, 3F12.3, F15.4)
    End Subroutine
