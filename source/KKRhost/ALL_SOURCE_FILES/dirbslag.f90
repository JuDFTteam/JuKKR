subroutine dirbslag(xi, y1i, y2i, y3i, y4i, y1, y2, y3, y4, ind1, n, imax)
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
  implicit none

! Dummy arguments
  integer :: imax, ind1, n
  real *8 :: xi
  real *8 :: y1i, y2i, y3i, y4i
  real *8 :: y1(imax), y2(imax), y3(imax), y4(imax)

! Local variables
  real *8 :: d, p, xd
  integer :: i, ind, inl, inu, j

  ind = ind1
  if (abs(xi-dble(ind))<1.0d-12) then
    y1i = y1(ind)
    y2i = y2(ind)
    y3i = y3(ind)
    y4i = y4(ind)
    return
  end if
! ------------------------------------- shift IND for irregular solution
  if (xi>dble(ind)) ind = ind + 1

  inl = max(1, ind-(n+1)/2)
  inu = inl + n - 1

  if (inu>imax) then
    inl = imax - n + 1
    inu = imax
  end if

  y1i = 0.0d0
  y2i = 0.0d0
  y3i = 0.0d0
  y4i = 0.0d0
  p = 1.0d0
  do j = inl, inu
    p = p*(xi-dble(j))
    d = 1.0d0
    do i = inl, inu
      if (i/=j) then
        xd = dble(j)
      else
        xd = xi
      end if
      d = d*(xd-dble(i))
    end do

    y1i = y1i + y1(j)/d
    y2i = y2i + y2(j)/d
    y3i = y3i + y3(j)/d
    y4i = y4i + y4(j)/d

  end do
  y1i = y1i*p
  y2i = y2i*p
  y3i = y3i*p
  y4i = y4i*p
end subroutine
