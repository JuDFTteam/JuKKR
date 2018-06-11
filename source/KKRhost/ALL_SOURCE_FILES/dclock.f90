! **********************************************************************
    Function dclock()
      Use mod_datatypes, Only: dp
      Real (Kind=dp) :: dclock
! **********************************************************************
!     .. External Functions ..
      Real (Kind=dp) :: etime, tarry(2)
      External :: etime

!     ..
      dclock = real(etime(tarry), kind=dp)
      Return
    End Function
