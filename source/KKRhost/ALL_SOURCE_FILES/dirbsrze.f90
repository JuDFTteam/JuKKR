    Subroutine dirbsrze(iest, xest, yest, yz, dy, nv, nuse)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *   diagonal rational function extrapolation to support the        *
!   *   Burlisch-Stoer method                                          *
!   *                                                                  *
!   *   see: numerical recipes chapter 15.4                            *
!   *                                                                  *
!   ********************************************************************

      Implicit None

! PARAMETER definitions
      Integer :: ncfmax, iseqmax, nusemax
      Parameter (ncfmax=8, iseqmax=30, nusemax=7)

! Dummy arguments
      Integer :: iest, nuse, nv
      Real (Kind=dp) :: xest
      Complex (Kind=dp) :: dy(nv), yest(nv), yz(nv)

! Local variables
      Complex (Kind=dp) :: b, b1, c, d(ncfmax, nusemax), ddy, v, yy
      Real (Kind=dp) :: fx(nusemax), x(iseqmax)
      Integer :: j, k, m1
      Save :: b, b1, c, d, ddy, fx, j, k, m1, v, x, yy

      x(iest) = xest
      If (iest==1) Then
        Do j = 1, nv
          yz(j) = yest(j)
          d(j, 1) = yest(j)
          dy(j) = yest(j)
        End Do
      Else
        m1 = min(iest, nuse)
        Do k = 1, m1 - 1
          fx(k+1) = x(iest-k)/xest
        End Do
        Do j = 1, nv
          yy = yest(j)
          v = d(j, 1)
          c = yy
          d(j, 1) = yy
          Do k = 2, m1
            b1 = fx(k)*v
            b = b1 - c
            If (b/=0._dp) Then
              b = (c-v)/b
              ddy = c*b
              c = b1*b
            Else
              ddy = v
            End If
            v = d(j, k)
            d(j, k) = ddy
            yy = yy + ddy
          End Do
          dy(j) = ddy
          yz(j) = yy
        End Do
      End If
    End Subroutine
