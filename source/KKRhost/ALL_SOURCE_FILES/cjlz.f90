    Function cjlz(l, z)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *   SPHERICAL BESSEL-FUNCTION  J(L,Z)  FOR COMPLEX ARGUMENT  Z     *
!   *                  see:  e.g. MERZBACHER EQ. (10.22)               *
!   *                                                                  *
!   ********************************************************************

      Implicit None

! PARAMETER definitions
      Complex (Kind=dp) :: c1
      Parameter (c1=(1.0E0_dp,0.0E0_dp))
      Integer :: lp2max
      Parameter (lp2max=25)


! Dummy arguments
      Integer :: l
      Complex (Kind=dp) :: z
      Complex (Kind=dp) :: cjlz

! Local variables
      Real (Kind=dp) :: dfac
      Complex (Kind=dp) :: dt, s(lp2max), t, zsq
      Integer :: i, k, llp1

      zsq = z*z
      llp1 = l + l + 1

      If (abs(zsq/real(llp1,kind=dp))<=10.E0_dp) Then

        dfac = 1.0E0_dp
        Do k = 3, llp1, 2
          dfac = dfac*real(k, kind=dp)
        End Do

        dt = c1
        t = c1
        Do i = 2, 400, 2
          dt = -dt*zsq/real(i*(i+llp1), kind=dp)
          t = t + dt
          If (abs(dt)<1.0E-10_dp) Go To 100
        End Do

100     Continue
        cjlz = t*z**l/dfac

      Else
        If (l>23) Stop '<cjlz>: l too large'

        s(2) = sin(z)/z
        If (l<=0) Then
          cjlz = s(2)*z**l
          Return
        End If

        s(1) = cos(z)
        Do i = 3, l + 2
          s(i) = (s(i-1)*(2*i-5)-s(i-2))/zsq
        End Do
        cjlz = s(l+2)*z**l

      End If
    End Function
