module mod_cint4pts

contains

subroutine cint4pts(y, jtop, z)
  ! ********************************************************************
  ! *                                                                  *
  ! *      perform the integral  Z(i)   =  INT   Y(i') di'             *
  ! *                                    R=0..R(i)                     *
  ! *                                                                  *
  ! *      via a 4-point integration formula                           *
  ! *                                                                  *
  ! *      JTOP:     Y is tabulated form 1 .. JTOP                     *
  ! *      Y(i):     function to be integrated                         *
  ! *                                                                  *
  ! *                       COMPLEX - VERSION                          *
  ! *                                                                  *
  ! ********************************************************************
  use :: mod_datatypes
  implicit none

  ! Dummy arguments
  integer :: jtop
  complex (kind=dp) :: y(jtop), z(jtop)

  ! Local variables
  integer :: i, ig, j, k, m, n1, n2
  real (kind=dp) :: q(5, 5), q5(5, 5)
  complex (kind=dp) :: s, svn

  data q5/0.e0_dp, 251.e0_dp, 232.e0_dp, 243.e0_dp, 224.e0_dp, 0.e0_dp, &
    646.e0_dp, 992.e0_dp, 918.e0_dp, 1024.e0_dp, 0.e0_dp, -264.e0_dp, &
    192.e0_dp, 648.e0_dp, 384.e0_dp, 0.e0_dp, 106.e0_dp, 32.e0_dp, 378.e0_dp, &
    1024.e0_dp, 0.e0_dp, -19.e0_dp, -8.e0_dp, -27.e0_dp, 224.e0_dp/

  do i = 1, 5
    do j = 1, 5
      q(i, j) = q5(i, j)/720.0e0_dp
    end do
  end do

  z(1) = cmplx(0.e0_dp, 0.e0_dp, kind=dp)
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

end subroutine cint4pts

end module mod_cint4pts
