subroutine cint4pts(y, jtop, z)
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
!   *                       COMPLEX - VERSION                          *
!   *                                                                  *
!   ********************************************************************
  use mod_DataTypes
  implicit none

! Dummy arguments
  integer :: jtop
  complex *16 :: y(jtop), z(jtop)

! Local variables
  integer :: i, ig, j, k, m, n1, n2
  real *8 :: q(5, 5), q5(5, 5)
  complex *16 :: s, svn

  data q5/0.d0, 251.d0, 232.d0, 243.d0, 224.d0, 0.d0, 646.d0, 992.d0, 918.d0, &
    1024.d0, 0.d0, -264.d0, 192.d0, 648.d0, 384.d0, 0.d0, 106.d0, 32.d0, &
    378.d0, 1024.d0, 0.d0, -19.d0, -8.d0, -27.d0, 224.d0/

  do i = 1, 5
    do j = 1, 5
      q(i, j) = q5(i, j)/720.0d0
    end do
  end do

  z(1) = cmplx(0.d0, 0.d0, kind=dp)
  svn = z(1)

  do ig = 1, jtop - 4, 4
    n1 = ig
    n2 = ig + 4
    do m = n1 + 1, n2
      i = m - n1 + 1
      s = svn
      do k = n1, n2
        j = k - n1 + 1
        s = s + q(i, j)*y(k)
      end do
      z(m) = s
    end do
    svn = z(n2)
  end do

  if (n2/=jtop) then
    n1 = jtop - 4
    n2 = jtop
    svn = z(n1)
    do m = n1 + 1, n2
      i = m - n1 + 1
      s = svn
      do k = n1, n2
        j = k - n1 + 1
        s = s + q(i, j)*y(k)
      end do
      z(m) = s
    end do
  end if

end subroutine
