    Subroutine potcut(imt1, irc1, ins, lmpot, r, vm2z, vspsme, vins, z1, irmd, &
      irmind)
      Use mod_datatypes, Only: dp
! **********************************************************************
! * set potential equal zero between muffin-tin and outer sphere       *
! **********************************************************************
      Implicit None
!     ..
!     .. Scalar Arguments ..
      Real (Kind=dp) :: z1
      Integer :: irmd, irmind
      Integer :: imt1, ins, irc1, lmpot
!     ..
!     .. Array Arguments ..
      Real (Kind=dp) :: r(*), vins(irmind:irmd, *), vm2z(*), vspsme(*)
!     ..
!     .. Local Scalars ..
      Integer :: ir, ist, lm
!     ..
!     .. Intrinsic Functions ..
      Intrinsic :: max
!     ..
      Write (1337, *) 'potcut: potential equal 2*Z/R between MT ', &
        'and outer sphere'
      Do ir = imt1 + 1, irc1
        vm2z(ir) = 2.0E0_dp*z1/r(ir)
        vspsme(ir) = 2.0E0_dp*z1/r(ir)
      End Do

      If (ins>=1) Then
        ist = max(irmind, imt1+1)
        Do ir = ist, irc1
          Do lm = 2, lmpot
            vins(ir, lm) = 0.0E0_dp
          End Do
        End Do
      End If
    End Subroutine ! SUBROUTINE POTCUT
