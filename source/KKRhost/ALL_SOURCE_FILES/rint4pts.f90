!------------------------------------------------------------------------------------
!> Summary: Perform an integral via a 4-point integration formula
!> Author: 
!> Perform an integral via a 4-point integration formula.
!> \begin{equation}
!> Z\left(i\right)=\int_{R=0}^{R\left(i\right)} Y\left(i'\right)di'
!> \end{equation}
!------------------------------------------------------------------------------------
module mod_rint4pts
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Perform an integral via a 4-point integration formula
  !> Author: 
  !> Category: numerical-tools, dirac, KKRhost
  !> Deprecated: False 
  !> Perform an integral via a 4-point integration formula.
  !> \begin{equation}
  !> Z\left(i\right)=\int_{R=0}^{R\left(i\right)} Y\left(i'\right)di'
  !> \end{equation}
  !-------------------------------------------------------------------------------
  subroutine rint4pts(y, jtop, z)

    implicit none

    ! Dummy arguments
    integer, intent(in) :: jtop !! Y is tabulated from 1 to JTOP
    real (kind=dp), dimension(jtop), intent(in) :: y  !! Function to be integrated
    real (kind=dp), dimension(jtop), intent(out) :: z !! Integrated function
    ! Local variables
    integer :: i, ig, j, k, m, n1, n2
    real (kind=dp) :: q(5, 5), q5(5, 5), s, svn
    data q5/0.e0_dp, 251.e0_dp, 232.e0_dp, 243.e0_dp, 224.e0_dp, 0.e0_dp, 646.e0_dp, 992.e0_dp, 918.e0_dp, 1024.e0_dp, 0.e0_dp, -264.e0_dp, 192.e0_dp, 648.e0_dp, 384.e0_dp, &
      0.e0_dp, 106.e0_dp, 32.e0_dp, 378.e0_dp, 1024.e0_dp, 0.e0_dp, -19.e0_dp, -8.e0_dp, -27.e0_dp, 224.e0_dp/

    do i = 1, 5
      do j = 1, 5
        q(i, j) = q5(i, j)/720.0e0_dp
      end do
    end do

    z(1) = 0.0e0_dp
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

  end subroutine rint4pts

end module mod_rint4pts
