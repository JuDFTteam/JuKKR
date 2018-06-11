    Subroutine csimpk(cf, cfint, ipan, ircut, drdi)
      Use mod_datatypes, Only: dp
!-----------------------------------------------------------------------
!     this subroutine does an integration up to rcut of an
!     complex function cf with an extended 3-point-simpson :

!                             rcut
!                      cfint = { cf(r') dr'
!                              0

!     modified for functions with kinks - at each kink the
!     integration is restarted .

!     attention : input cf is destroyed !

!-----------------------------------------------------------------------
!.. Scalar Arguments ..
      Complex (Kind=dp) :: cfint
      Integer :: ipan
!..
!.. Array Arguments ..
      Complex (Kind=dp) :: cf(*)
      Real (Kind=dp) :: drdi(*)
      Integer :: ircut(0:ipan)
!..
!.. Local Scalars ..
      Real (Kind=dp) :: a1, a2
      Integer :: i, ien, ip, ist, n
!..
!.. External Functions ..
      Complex (Kind=dp) :: csum
      External :: csum
!..
!.. Intrinsic Functions ..
      Intrinsic :: mod
!     ..
      a1 = 4.0E0_dp/3.0E0_dp
      a2 = 2.0E0_dp/3.0E0_dp
      cfint = 0.0E0_dp

      Do ip = 1, ipan

!---> loop over kinks

        ist = ircut(ip-1) + 1
        ien = ircut(ip)

        Do i = ist, ien
          cf(i) = cf(i)*drdi(i)
        End Do

        If (mod(ien-ist,2)==0) Then
          cfint = cfint + (cf(ist)-cf(ien))/3.0E0_dp
          ist = ist + 1
          n = (ien-ist+1)/2

        Else
!---> four point lagrange integration for the first step
          cfint = cfint + (9.0E0_dp*cf(ist)+19.0E0_dp*cf(ist+1)-5.0E0_dp*cf( &
            ist+2)+cf(ist+3))/24.0E0_dp + (cf(ist+1)-cf(ien))/3.0E0_dp
          ist = ist + 2
          n = (ien-ist+1)/2
        End If

!---> calculate with an extended 3-point-simpson

        cfint = cfint + a1*csum(n, cf(ist), 2) + a2*csum(n, cf(ist+1), 2)
      End Do

    End Subroutine
