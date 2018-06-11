    Subroutine force(flm, flmc, lmax, nspin, nstart, nend, rhoc, v, r, drdi, &
      irws)
      Use mod_datatypes, Only: dp
!-----------------------------------------------------------------------
!     calculates the force on nucleus m
!     from a given non spherical charge density at the nucleus site r
!     with core correction (coulomb contribution)

!-----------------------------------------------------------------------
      Implicit None

!     .. Parameters ..
      Include 'inc.p'
      Integer :: lmpotd
      Parameter (lmpotd=(lpotd+1)**2)
!..
!.. Local Scalars ..
      Integer :: lmax, nend, nspin, nstart
!..
!.. Local Arrays ..
      Real (Kind=dp) :: drdi(irmd, *), flm(-1:1, *), flmc(-1:1, *), &
        r(irmd, *), rhoc(irmd, *), v(irmd, lmpotd, *)
      Integer :: irws(*)
!..
!.. External Subroutines ..
      Real (Kind=dp) :: dv, fac, pi, rws, vint1
      Integer :: i, iatyp, ipot, irws1, ispin, lm, m
!..
!.. Save statement ..
      Real (Kind=dp) :: flmh(-1:1, natypd), v1(irmd)
!..
!.. Intrinsic Functions ..
      External :: simp3
!..

      Save :: pi

!---> loop over rep. atoms
      Intrinsic :: atan, sqrt

      pi = 4.E0_dp*atan(1.E0_dp)
      fac = sqrt((4.0E0_dp*pi)/3.0E0_dp)
      If (lmax<1) Then
        Write (6, Fmt=100)
        Stop

      End If



      Do iatyp = nstart, nend


        irws1 = irws(iatyp)
        rws = r(irws1, iatyp)
!---> initialize v1


        Do m = -1, 1
          lm = 2 + m + 1

!---> determine the right potential numbers

          Do i = 1, irws1
            v1(i) = 0.0E0_dp
          End Do

          Do ispin = 1, nspin
!---> determine the derivative of the potential using a 5-point formular


            ipot = nspin*(iatyp-1) + ispin



            dv = (-3.0E0_dp*v(1,lm,ipot)-10.0E0_dp*v(2,lm,ipot)+ &
              18.0E0_dp*v(3,lm,ipot)-6.0E0_dp*v(4,lm,ipot)+v(5,lm,ipot))/ &
              (12.0E0_dp*drdi(2,iatyp))

            v1(2) = rhoc(2, ipot)*(2.0E0_dp*v(2,lm,ipot)/r(2,iatyp)+dv)/ &
              (4.0E0_dp*pi) + v1(2)

            Do i = 3, irws1 - 2

              dv = (v(i-2,lm,ipot)-v(i+2,lm,ipot)+8.0E0_dp*(v(i+1,lm, &
                ipot)-v(i-1,lm,ipot)))/(12.0E0_dp*drdi(i,iatyp))

              v1(i) = rhoc(i, ipot)*(2.0E0_dp*v(i,lm,ipot)/r(i,iatyp)+dv)/ &
                (4.0E0_dp*pi) + v1(i)
            End Do
!---> integrate with simpson subroutine
            dv = (-v(irws1-4,lm,ipot)+6.0E0_dp*v(irws1-3,lm,ipot)- &
              18.0E0_dp*v(irws1-2,lm,ipot)+10.0E0_dp*v(irws1-1,lm,ipot)+ &
              3.0E0_dp*v(irws1,lm,ipot))/(12.0E0_dp*drdi(irws1-1,iatyp))
            v1(irws1-1) = rhoc(irws1-1, ipot)*(2.0E0_dp*v(irws1-1,lm,ipot)/r( &
              irws1-1,iatyp)+dv)/(4.0E0_dp*pi) + v1(irws1-1)

            dv = (3.0E0_dp*v(irws1-4,lm,ipot)-16.0E0_dp*v(irws1-3,lm,ipot)+ &
              36.0E0_dp*v(irws1-2,lm,ipot)-48.0E0_dp*v(irws1-1,lm,ipot)+ &
              25.0E0_dp*v(irws1,lm,ipot))/(12.0E0_dp*drdi(irws1,iatyp))

            v1(irws1) = rhoc(irws1, ipot)*(2.0E0_dp*v(irws1,lm,ipot)/r(irws1, &
              iatyp)+dv)/(4.0E0_dp*pi) + v1(irws1)
          End Do



          Call simp3(v1, vint1, 1, irws1, drdi(1,iatyp))

          flmh(m, iatyp) = fac*flm(m, iatyp)
          flmc(m, iatyp) = -fac*vint1
          flm(m, iatyp) = flmh(m, iatyp) + flmc(m, iatyp)


        End Do
!-----------------------------------------------------------------------
!     calculates the force on nucleus m
      End Do
!     from a given non spherical charge density at the nucleus site r
100   Format (13X, 'error stop in subroutine force :', &
        ' the charge density has to contain non spherical', &
        ' contributions up to l=1 at least ')
!     with core correction (coulomb contribution)
    End Subroutine
