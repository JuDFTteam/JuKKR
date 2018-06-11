    Subroutine csinwd(f, fint, lmmsqd, irmind, irmd, irmin, ipan, ircut)
      Use mod_datatypes, Only: dp
!-----------------------------------------------------------------------
!     this subroutine does an inwards integration of llmax
!     functions f with an extended 3-point-simpson :


!                               irmax
!                   fint(ll,i) = { f(ll,i') di'
!                                ir

!  the starting value for this integration at ist - 1 is determined by
!    a 4 point lagrangian integration  , coefficients given by
!    m. abramowitz and i.a. stegun, handbook of mathematical functions,
!    nbs applied mathematics series 55 (1968)

!  attention in case of radial integration :
!       the weights drdi have to be multiplied before calling this
!       subroutine .

!                                     b. drittler mar. 1989

!    modified for functions with kinks - at each kink the integration
!      is restarted

!    attention : it is supposed that irmin + 3 is less than imt !


!                                     b. drittler july 1989
!    modified by m. ogura, june 2015
!-----------------------------------------------------------------------
!     ..
!.. Parameters ..
      Real (Kind=dp) :: a1, a2, a3
      Parameter (a1=5.E0_dp/12.E0_dp, a2=8.E0_dp/12.E0_dp, &
        a3=-1.E0_dp/12.E0_dp)
!..
!.. Scalar Arguments ..
      Integer :: ipan, irmd, irmind, lmmsqd, irmin
!..
!.. Array Arguments ..
      Complex (Kind=dp) :: f(lmmsqd, irmind:irmd), fint(lmmsqd, irmind:irmd)
      Integer :: ircut(0:ipan)
!..
!.. Local Scalars ..
      Integer :: i, ien, ip, ist, ll

!---> loop over kinks

      Do ip = ipan, 1, -1
        ist = ircut(ip)
        ien = ircut(ip-1) + 1
        If (ip==1) ien = irmin

        If (ip==ipan) Then
          Do ll = 1, lmmsqd
            fint(ll, ist) = 0.0E0_dp
          End Do

        Else
          Do ll = 1, lmmsqd
            fint(ll, ist) = fint(ll, ist+1)
          End Do
        End If

!---> calculate fint with an extended 3-point-simpson

        Do i = ist, ien + 2, -2
          Do ll = 1, lmmsqd
            fint(ll, i-1) = fint(ll, i) + f(ll, i)*a1 + f(ll, i-1)*a2 + &
              f(ll, i-2)*a3
            fint(ll, i-2) = fint(ll, i-1) + f(ll, i)*a3 + f(ll, i-1)*a2 + &
              f(ll, i-2)*a1
          End Do
        End Do
        If (mod(ist-ien,2)==1) Then
          Do ll = 1, lmmsqd
            fint(ll, ien) = fint(ll, ien+1) + f(ll, ien)*a1 + &
              f(ll, ien+1)*a2 + f(ll, ien+2)*a3
          End Do
        End If
      End Do

    End Subroutine
