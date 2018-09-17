module mod_spatpr
  use :: mod_datatypes, only: dp
  private :: dp

contains

  ! ************************************************************************
  subroutine spatpr(a, b, c, v)
    ! ************************************************************************
    ! SPATPR COMPUTES THE SPATIAL PRODUCT OF THREE VECTORS A,B AND C
    ! RETURNING IT INTO V: V=AXB.C.
    ! ------------------------------------------------------------------------


    real (kind=dp), intent (in) :: a(*)
    real (kind=dp), intent (in) :: b(*)
    real (kind=dp), intent (in) :: c(*)
    real (kind=dp), intent (out) :: v

    v = 0.0e0_dp
    v = v + c(1)*(a(2)*b(3)-a(3)*b(2))
    v = v + c(2)*(a(3)*b(1)-a(1)*b(3))
    v = v + c(3)*(a(1)*b(2)-a(2)*b(1))
    return
  end subroutine spatpr

end module mod_spatpr
