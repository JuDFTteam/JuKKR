    Subroutine simp3(f, fint, istart, iend, drdi)
      Use mod_datatypes, Only: dp
!-----------------------------------------------------------------------
!     this subroutine does an integration from istart to iend of
!     the real function f with an extended 3-point-simpson :

!                          r(istart)

!                       fint = { f(r') dr'

!                           r(iend)

!-----------------------------------------------------------------------
!.. Scalar Arguments ..
      Real (Kind=dp) :: fint
      Integer :: iend, istart
!..
!.. Array Arguments ..
      Real (Kind=dp) :: drdi(*), f(*)
!..
!.. Local Scalars ..
      Real (Kind=dp) :: a1, a2
      Integer :: i, ist
!..
!.. Intrinsic Functions ..
      Intrinsic :: mod
!..
      a1 = 4.0E0_dp/3.0E0_dp
      a2 = 2.0E0_dp/3.0E0_dp

!---> initialize fint

      If (mod(iend-istart,2)==0) Then
        fint = f(istart)*drdi(istart)/3.0E0_dp
        ist = istart + 1

      Else
        fint = (f(istart+3)*drdi(istart+3)-5.0E0_dp*f(istart+2)*drdi(istart+2) &
          +19.0E0_dp*f(istart+1)*drdi(istart+1)+9.0E0_dp*f(istart)*drdi(istart &
          ))/24.0E0_dp + f(istart+1)*drdi(istart+1)/3.0E0_dp
        ist = istart + 2
      End If

!---> calculate with an extended 3-point-simpson

      Do i = ist, iend - 1, 2
        fint = fint + a1*f(i)*drdi(i) + a2*f(i+1)*drdi(i+1)
      End Do
      fint = fint - f(iend)*drdi(iend)/3.0E0_dp

    End Subroutine
