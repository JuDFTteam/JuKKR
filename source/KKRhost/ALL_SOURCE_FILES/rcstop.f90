! ************************************************************************
    Subroutine rcstop(c)
! ************************************************************************
!     .. Scalar Arguments ..

      Character (Len=8), Intent (In) :: c

!     ..
      Print *, 'ERROR: STOP AT POSITION ', c
      Stop

    End Subroutine
