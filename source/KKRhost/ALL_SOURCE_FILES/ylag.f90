    Function ylag(xi, x, y, ind1, n1, imax)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   * lagrangian interpolation                                         *
!   * xi is interpolated entry into x-array                            *
!   * n is the order of lagrangran interpolation                       *
!   * y is array from which ylag is obtained by interpolation          *
!   * ind is the min-i for x(i).gt.xi                                  *
!   * if ind=0,x-array will be searched                                *
!   * imax is max index of x-and y-arrays                              *
!   *                                                                  *
!   * 07/12/94  HE  arg. IEX removed                                   *
!   ********************************************************************

      Implicit None

! Dummy arguments
      Integer :: imax, ind1, n1
      Real (Kind=dp) :: xi
      Real (Kind=dp) :: x(imax), y(imax)
      Real (Kind=dp) :: ylag

! Local variables
      Real (Kind=dp) :: d, p, s, xd
      Integer :: i, ind, inl, inu, j, n
      Save :: d, i, ind, inl, inu, j, n, p, s, xd

      ind = ind1
      n = n1
      If (n>imax) n = imax
      If (ind>0) Go To 110
      Do j = 1, imax
        If (abs(xi-x(j))<1.0E-12_dp) Go To 150
        If (xi<x(j)) Go To 100
        If (xi==x(j)) Go To 150
      End Do
      Go To 120
100   Continue
      ind = j
110   Continue
      If (ind>1) Then
      End If
      inl = ind - (n+1)/2
      If (inl<=0) inl = 1
      inu = inl + n - 1
      If (inu<=imax) Go To 130
120   Continue
      inl = imax - n + 1
      inu = imax
130   Continue
      s = 0.0E0_dp
      p = 1.0E0_dp
      Do j = inl, inu
        p = p*(xi-x(j))
        d = 1.0E0_dp
        Do i = inl, inu
          If (i/=j) Then
            xd = x(j)
          Else
            xd = xi
          End If
          d = d*(xd-x(i))
        End Do
        s = s + y(j)/d
      End Do
      ylag = s*p
140   Continue
      Return
150   Continue
      ylag = y(j)
      Go To 140
    End Function
