    Subroutine rint4pts(y, jtop, z)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *      perform the integral  Z(i)   =  INT   Y(i') di'             *
!   *                                    R=0..R(i)                     *
!   *                                                                  *
!   *      via a 4-point integration formula                           *
!   *                                                                  *
!   *      JTOP:     Y is tabulated form 1 .. JTOP                     *
!   *      Y(i):     function to be integrated                         *
!   *                                                                  *
!   *                       REAL    - VERSION                          *
!   *                                                                  *
!   ********************************************************************
      Implicit None

! Dummy arguments
      Integer :: jtop
      Real (Kind=dp) :: y(jtop), z(jtop)
! Local variables
      Integer :: i, ig, j, k, m, n1, n2
      Real (Kind=dp) :: q(5, 5), q5(5, 5), s, svn
      Data q5/0.E0_dp, 251.E0_dp, 232.E0_dp, 243.E0_dp, 224.E0_dp, 0.E0_dp, &
        646.E0_dp, 992.E0_dp, 918.E0_dp, 1024.E0_dp, 0.E0_dp, -264.E0_dp, &
        192.E0_dp, 648.E0_dp, 384.E0_dp, 0.E0_dp, 106.E0_dp, 32.E0_dp, &
        378.E0_dp, 1024.E0_dp, 0.E0_dp, -19.E0_dp, -8.E0_dp, -27.E0_dp, &
        224.E0_dp/

      Do i = 1, 5
        Do j = 1, 5
          q(i, j) = q5(i, j)/720.0E0_dp
        End Do
      End Do

      z(1) = 0.0E0_dp
      svn = z(1)

      Do ig = 1, jtop - 4, 4
        n1 = ig
        n2 = ig + 4
        Do m = n1 + 1, n2
          i = m - n1 + 1
          s = svn
          Do k = n1, n2
            j = k - n1 + 1
            s = s + q(i, j)*y(k)
          End Do
          z(m) = s
        End Do
        svn = z(n2)
      End Do

      If (n2/=jtop) Then
        n1 = jtop - 4
        n2 = jtop
        svn = z(n1)
        Do m = n1 + 1, n2
          i = m - n1 + 1
          s = svn
          Do k = n1, n2
            j = k - n1 + 1
            s = s + q(i, j)*y(k)
          End Do
          z(m) = s
        End Do
      End If

    End Subroutine
