subroutine corlsd(rs, zta, ec, vcup, vcdn, ecrs, eczta, alfc)
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
  double precision :: alfc, ec, ecrs, eczta, rs, vcdn, vcup, zta
!..
!.. Local Scalars ..
  double precision :: alfm, alfrsm, comm, ep, eprs, eu, eurs, f, fz, fzz, gam, &
    thrd, thrd4, z4
!..
!.. External Subroutines ..
  external :: gcor91
!..
!.. Save statement ..
  save :: gam, fzz, thrd, thrd4
!..
!.. Data statements ..
!.....-----------------------------------------------------------------
  data gam, fzz/0.5198421d0, 1.709921d0/
  data thrd, thrd4/0.333333333333d0, 1.333333333333d0/
!..
!.....-----------------------------------------------------------------
  f = ((1.d0+zta)**thrd4+(1.d0-zta)**thrd4-2.d0)/gam
  call gcor91(0.0310907d0, 0.21370d0, 7.5957d0, 3.5876d0, 1.6382d0, 0.49294d0, &
    1.00d0, rs, eu, eurs)
  call gcor91(0.01554535d0, 0.20548d0, 14.1189d0, 6.1977d0, 3.3662d0, &
    0.62517d0, 1.00d0, rs, ep, eprs)
  call gcor91(0.0168869d0, 0.11125d0, 10.357d0, 3.6231d0, 0.88026d0, &
    0.49671d0, 1.00d0, rs, alfm, alfrsm)
!  alfm is minus the spin stiffness alfc
  alfc = -alfm
  z4 = zta**4
  ec = eu*(1.d0-f*z4) + ep*f*z4 - alfm*f*(1.d0-z4)/fzz
!  energy done. now the potential:
  ecrs = eurs*(1.d0-f*z4) + eprs*f*z4 - alfrsm*f*(1.d0-z4)/fzz
  fz = thrd4*((1.d0+zta)**thrd-(1.d0-zta)**thrd)/gam
  eczta = 4.d0*(zta**3)*f*(ep-eu+alfm/fzz) + fz*(z4*ep-z4*eu-(1.d0-z4)*alfm/ &
    fzz)
  comm = ec - rs*ecrs/3.d0 - zta*eczta
  vcup = comm + eczta
  vcdn = comm - eczta

  return
end subroutine
