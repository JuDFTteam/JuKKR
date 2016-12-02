!module sphj_mod
!implicit none
!private
!public :: sphj
!contains


! from the specfun package as included in scipy
  subroutine sphj(n,x,nm,sj,dj)
!       modified to allow n=0 case (also in csphjy, sphy)
!       =======================================================
!       purpose: compute spherical bessel functions jn(x) and
!                their derivatives
!       input :  x --- argument of jn(x)
!                n --- order of jn(x)  ( n = 0,1,... )
!       output:  sj(n) --- jn(x)
!                dj(n) --- jn'(x)
!                nm --- highest order computed
!       routines called:
!                msta1 and msta2 for computing the starting
!                point for backward recurrence
!       =======================================================
    implicit none
    integer, intent(in) :: n
    integer, intent(out) :: nm
    double precision, intent(in) :: x
    double precision, intent(out) :: sj(0:n), dj(0:n)
    
    integer, external :: msta1, msta2
    double precision :: f, f0, f1, cs, sa, sb
    integer :: k, m
    nm = n
    if (dabs(x) < 1.d-100) then
      do k = 0, n
        sj(k) = 0.d0
        dj(k) = 0.d0
      enddo ! k
      sj(0) = 1.0d0
      if (n > 0) dj(1) = .3333333333333333d0
      return
    endif
    sj(0) = dsin(x)/x
    dj(0) = (dcos(x) - dsin(x)/x)/x
    if (n < 1) return
    sj(1) = (sj(0)-dcos(x))/x
    if (n >= 2) then
      sa = sj(0)
      sb = sj(1)
      m = msta1(x, 200)
      if (m < n) then
        nm = m
      else
        m = msta2(x, n, 15)
      endif
      f = 0.d0
      f0 = 0.d0
      f1 = 1.d0 - 100
      do k = m, 0, -1
        f = (2*k + 3.d0)*f1/x - f0
        if (k <= nm) sj(k) = f
        f0 = f1
        f1 = f
      enddo ! k
      cs = 0.d0
      if (dabs(sa) > dabs(sb)) then
        cs = sa/f
      else
        cs = sb/f0
      endif
      do k= 0, nm
        sj(k) = cs*sj(k)
      enddo ! k
    endif
    do k = 1, nm
      dj(k) = sj(k - 1) - (k + 1.d0)*sj(k)/x
    enddo ! k
  endsubroutine ! sphj

  integer function msta1(x,mp)
!       ===================================================
!       purpose: determine the starting point for backward
!                recurrence such that the magnitude of
!                jn(x) at that point is about 10^(-mp)
!       input :  x     --- argument of jn(x)
!                mp    --- value of magnitude
!       output:  msta1 --- starting point
!       ===================================================
    implicit none
    double precision, intent(in) :: x
    integer, intent(in) :: mp
    double precision, external :: envj
    double precision :: a0, f0, f1, f
    integer :: n0, n1, nn, it
    
    a0 = dabs(x)
    n0 = int(1.1d0*a0) + 1
    f0 = envj(n0,a0) - mp
    n1 = n0 + 5
    f1 = envj(n1, a0) - mp
    do it = 1, 20
      nn = int(n1 - (n1 - n0)/(1.d0 - f0/f1))
      f = envj(nn, a0) - mp
      if (abs(nn - n1) < 1) exit
      n0 = n1
      f0 = f1
      n1 = nn
      f1 = f
    enddo ! it
    msta1 = nn
  endfunction ! msta1

  integer function msta2(x,n,mp)
!       ===================================================
!       purpose: determine the starting point for backward
!                recurrence such that all jn(x) has mp
!                significant digits
!       input :  x  --- argument of jn(x)
!                n  --- order of jn(x)
!                mp --- significant digit
!       output:  msta2 --- starting point
!       ===================================================
    implicit none
    double precision, intent(in) :: x
    integer, intent(in) :: n, mp
    double precision, external :: envj
    double precision :: a0, hmp, ejn, obj, f0, f1, f
    integer :: n0, n1, nn, it
    a0 = dabs(x)
    hmp = 0.5d0*mp
    ejn = envj(n, a0)
    if (ejn <= hmp) then
      obj = mp
      n0 = int(1.1*a0) + 1
    else
      obj = hmp + ejn
      n0 = n
    endif
    f0 = envj(n0, a0) - obj
    n1 = n0 + 5
    f1 = envj(n1, a0) - obj
    do it = 1, 20
      nn = int(n1 - (n1 - n0)/(1.d0 - f0/f1))
      f = envj(nn, a0) - obj
      if (abs(nn - n1) < 1) exit
      n0 = n1
      f0 = f1
      n1 = nn
      f1 = f
    enddo ! it
    msta2 = nn + 10
  endfunction ! msta2

  double precision function envj(n,x)
  implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x
    envj = 0.5d0*dlog10(6.28d0*n) - n*dlog10(1.36d0*x/n)
  endfunction ! envj

!endmodule