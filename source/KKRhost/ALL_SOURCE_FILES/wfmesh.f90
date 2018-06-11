    Subroutine wfmesh(e, ek, cvlight, nsra, z, r, s, rs, irm, irmd, lmaxd)
      Use mod_datatypes, Only: dp
!.. Scalar Arguments ..
      Complex (Kind=dp) :: e, ek
      Real (Kind=dp) :: cvlight, z
      Integer :: irm, irmd, lmaxd, nsra
!..
!.. Intrinsic Functions ..
      Intrinsic :: real, sqrt
!..
!.. Array Arguments ..
      Real (Kind=dp) :: r(irmd), rs(irmd, 0:lmaxd), s(0:lmaxd)
!..
!.. Local Scalars ..
      Real (Kind=dp) :: s1
      Integer :: ir, l
!..
      If (nsra==1) ek = sqrt(e)
      If (nsra==2) ek = sqrt(e+e*e/(cvlight*cvlight))
      Do l = 0, lmaxd

        If (nsra==2) Then
          s1 = sqrt(real(l*l+l+1,kind=dp)-4.0E0_dp*z*z/(cvlight*cvlight))
          If (z==0.0E0_dp) s1 = real(l, kind=dp)
        Else
          s1 = real(l, kind=dp)
        End If
        s(l) = s1
        rs(1, l) = 0.0E0_dp
        Do ir = 2, irm
          rs(ir, l) = r(ir)**s1
        End Do
        Do ir = irm + 1, irmd
          rs(ir, l) = 0.0E0_dp
        End Do

      End Do
      Return
    End Subroutine
