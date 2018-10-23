module mod_gcor91
  
  private
  public :: gcor91

contains

  !-------------------------------------------------------------------------------
  !> Summary: Helper for PW91 correlation
  !> Author: 
  !> Category: KKRhost, xc-potential
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Called by corlsd
  !-------------------------------------------------------------------------------
  subroutine gcor91(a, a1, b1, b2, b3, b4, p, rs, gg, ggrs)

    use :: mod_datatypes, only: dp
    implicit none

    real (kind=dp) :: a, a1, b1, b2, b3, b4, gg, ggrs, p, rs
    real (kind=dp) :: p1, q0, q1, q2, q3, rs12, rs32, rsp

    intrinsic :: log, sqrt


    p1 = p + 1.e0_dp
    q0 = -2.e0_dp*a*(1.e0_dp+a1*rs)
    rs12 = sqrt(rs)
    rs32 = rs12**3
    rsp = rs**p
    q1 = 2.e0_dp*a*(b1*rs12+b2*rs+b3*rs32+b4*rs*rsp)
    q2 = log(1.e0_dp+1.e0_dp/q1)
    gg = q0*q2
    q3 = a*(b1/rs12+2.e0_dp*b2+3.e0_dp*b3*rs12+2.e0_dp*b4*p1*rsp)
    ggrs = -2.e0_dp*a*a1*q2 - q0*q3/(q1**2+q1)
    return

  end subroutine gcor91

end module mod_gcor91
