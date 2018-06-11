    Subroutine corlsd(rs, zta, ec, vcup, vcdn, ecrs, eczta, alfc)
      Use mod_datatypes, Only: dp
!.....-----------------------------------------------------------------
!     uniform-gas correlation of perdew and wang 1991
!.....-----------------------------------------------------------------
!     input: seitz radius (rs), relative spin polarization (zta)
!     output: correlation energy per electron (ec),
!             up- and down-spin potentials (vcup,vcdn),
!             derivatives of ec wrt rs (ecrs) &zta (eczta).
!     output: correlation contribution (alfc) to the spin stiffness
!.....-----------------------------------------------------------------
!.. Scalar Arguments ..
      Real (Kind=dp) :: alfc, ec, ecrs, eczta, rs, vcdn, vcup, zta
!..
!.. Local Scalars ..
      Real (Kind=dp) :: alfm, alfrsm, comm, ep, eprs, eu, eurs, f, fz, fzz, &
        gam, thrd, thrd4, z4
!..
!.. External Subroutines ..
      External :: gcor91
!..
!.. Save statement ..
      Save :: gam, fzz, thrd, thrd4
!..
!.. Data statements ..
!.....-----------------------------------------------------------------
      Data gam, fzz/0.5198421E0_dp, 1.709921E0_dp/
      Data thrd, thrd4/0.333333333333E0_dp, 1.333333333333E0_dp/
!..
!.....-----------------------------------------------------------------
      f = ((1.E0_dp+zta)**thrd4+(1.E0_dp-zta)**thrd4-2.E0_dp)/gam
      Call gcor91(0.0310907E0_dp, 0.21370E0_dp, 7.5957E0_dp, 3.5876E0_dp, &
        1.6382E0_dp, 0.49294E0_dp, 1.00E0_dp, rs, eu, eurs)
      Call gcor91(0.01554535E0_dp, 0.20548E0_dp, 14.1189E0_dp, 6.1977E0_dp, &
        3.3662E0_dp, 0.62517E0_dp, 1.00E0_dp, rs, ep, eprs)
      Call gcor91(0.0168869E0_dp, 0.11125E0_dp, 10.357E0_dp, 3.6231E0_dp, &
        0.88026E0_dp, 0.49671E0_dp, 1.00E0_dp, rs, alfm, alfrsm)
!  alfm is minus the spin stiffness alfc
      alfc = -alfm
      z4 = zta**4
      ec = eu*(1.E0_dp-f*z4) + ep*f*z4 - alfm*f*(1.E0_dp-z4)/fzz
!  energy done. now the potential:
      ecrs = eurs*(1.E0_dp-f*z4) + eprs*f*z4 - alfrsm*f*(1.E0_dp-z4)/fzz
      fz = thrd4*((1.E0_dp+zta)**thrd-(1.E0_dp-zta)**thrd)/gam
      eczta = 4.E0_dp*(zta**3)*f*(ep-eu+alfm/fzz) + fz*(z4*ep-z4*eu-(1.E0_dp- &
        z4)*alfm/fzz)
      comm = ec - rs*ecrs/3.E0_dp - zta*eczta
      vcup = comm + eczta
      vcdn = comm - eczta

      Return
    End Subroutine
