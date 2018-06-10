! **********************************************************************
double precision function dclock()
! **********************************************************************
!     .. External Functions ..
  real :: etime, tarry(2)
  external :: etime

!     ..
  dclock = dble(etime(tarry))
  return
end function
