    Subroutine forceh(cmom, flmh, lmax, nspin, nstart, nend, r2rho, v, r, &
      drdi, irws, z)
      Use mod_datatypes, Only: dp
!-----------------------------------------------------------------------
!     calculates the force on nucleus m with hellmann - feynman theorem
!     from a given non spherical charge density at the nucleus site r


!-----------------------------------------------------------------------
      Implicit None
!.. Parameters ..
      Include 'inc.p'
      Integer :: lmpotd
      Parameter (lmpotd=(lpotd+1)**2)
!..
!.. Local Scalars ..
      Integer :: lmax, nend, nspin, nstart
!..
!.. Local Arrays ..
      Real (Kind=dp) :: cmom(lmpotd, *), drdi(irmd, *), flmh(-1:1, *), &
        r(irmd, *), r2rho(irmd, lmpotd, natypd, *), v(irmd, lmpotd, *), z(*)
      Integer :: irws(*)
!..
!.. External Subroutines ..
      Real (Kind=dp) :: pi, rws, vint1
      Integer :: i, iatyp, ipot, irws1, lm, m
!..
!.. Save statement ..
      Real (Kind=dp) :: flm(-1:1, 2), v1(irmd)
!..

      External :: simp3
!.. Intrinsic Functions ..
!..
      Save :: pi


!---> loop over the rep. atoms
      Intrinsic :: atan

      pi = 4.E0_dp*atan(1.E0_dp)
      If (lmax<1) Then
        Write (6, Fmt=100)
        Stop

      End If
!---> reading the right Wigner-S. radius


      Do iatyp = nstart, nend
!---> determine the right potential numbers


        irws1 = irws(iatyp)
        rws = r(irws1, iatyp)


!---> integrate with simpson subroutine
        ipot = nspin*(iatyp-1) + 1

        Do m = -1, 1
          lm = 2 + m + 1

          v1(1) = 0.0E0_dp
          Do i = 2, irws1
            v1(i) = r2rho(i, lm, iatyp, 1)*(r(i,iatyp)**(-2.0E0_dp))
          End Do

!---> use coulomb potential to determine extra atomic contribution

          Call simp3(v1, vint1, 1, irws1, drdi(1,iatyp))

          flm(m, 1) = 2.0E0_dp*vint1
!---> total Hellman-Feynman force


          flm(m, 2) = v(irws1, lm, ipot)*(3.0E0_dp/(4.0E0_dp*pi*rws)) - &
            2.0E0_dp*cmom(lm, iatyp)/(rws**3)


!-----------------------------------------------------------------------
          flmh(m, iatyp) = (flm(m,1)+flm(m,2))*z(iatyp)
        End Do
      End Do
!     calculates the force on nucleus m with hellmann - feynman theorem
!     from a given non spherical charge density at the nucleus site r
100   Format (13X, 'error stop in subroutine force :', &
        ' the charge density has to contain non spherical', &
        ' contributions up to l=1 at least ')

    End Subroutine
