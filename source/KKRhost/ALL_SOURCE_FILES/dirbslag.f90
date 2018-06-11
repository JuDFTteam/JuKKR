    Subroutine dirbslag(xi, y1i, y2i, y3i, y4i, y1, y2, y3, y4, ind1, n, imax)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *      lagrangian interpolation of Y(X) at position XI             *
!   *                                                                  *
!   *      XI      entry into x-array                                  *
!   *              for regular solution:   X(IND1-1) < XI <=X(IND1)    *
!   *              for irregular solution: X(IND1)   < XI <=X(IND1+1)  *
!   *      X/Y     X/Y-arrays                                          *
!   *      N       order of lagrangian interpolation                   *
!   *      IND     min-I for which  X(I) > XI                          *
!   *      IMAX    max index of X/Y-arrays                             *
!   *                                                                  *
!   ********************************************************************
      Implicit None

! Dummy arguments
      Integer :: imax, ind1, n
      Real (Kind=dp) :: xi
      Real (Kind=dp) :: y1i, y2i, y3i, y4i
      Real (Kind=dp) :: y1(imax), y2(imax), y3(imax), y4(imax)

! Local variables
      Real (Kind=dp) :: d, p, xd
      Integer :: i, ind, inl, inu, j

      ind = ind1
      If (abs(xi-real(ind,kind=dp))<1.0E-12_dp) Then
        y1i = y1(ind)
        y2i = y2(ind)
        y3i = y3(ind)
        y4i = y4(ind)
        Return
      End If
! ------------------------------------- shift IND for irregular solution
      If (xi>real(ind,kind=dp)) ind = ind + 1

      inl = max(1, ind-(n+1)/2)
      inu = inl + n - 1

      If (inu>imax) Then
        inl = imax - n + 1
        inu = imax
      End If

      y1i = 0.0E0_dp
      y2i = 0.0E0_dp
      y3i = 0.0E0_dp
      y4i = 0.0E0_dp
      p = 1.0E0_dp
      Do j = inl, inu
        p = p*(xi-real(j,kind=dp))
        d = 1.0E0_dp
        Do i = inl, inu
          If (i/=j) Then
            xd = real(j, kind=dp)
          Else
            xd = xi
          End If
          d = d*(xd-real(i,kind=dp))
        End Do

        y1i = y1i + y1(j)/d
        y2i = y2i + y2(j)/d
        y3i = y3i + y3(j)/d
        y4i = y4i + y4(j)/d

      End Do
      y1i = y1i*p
      y2i = y2i*p
      y3i = y3i*p
      y4i = y4i*p
    End Subroutine
