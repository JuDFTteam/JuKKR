! ************************************************************************
      subroutine cinit(n, a)
! ************************************************************************
!-----------------------------------------------------------------------
!     initialize the first n values of a complex array a with zero
!-----------------------------------------------------------------------
!     .. scalar arguments ..
      integer, intent(in) :: n
!     ..
!     .. array arguments ..
      double complex, intent(out) :: a(*)
!     ..
!     .. local scalars ..
      integer :: i
      double complex, parameter :: zero=(0.d0, 0.d0)
!     ..
      do i = 1, n
        a(i) = zero
      enddo ! i
      endsubroutine ! cinit
