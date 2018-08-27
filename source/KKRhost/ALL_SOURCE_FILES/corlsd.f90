module mod_corlsd

contains

subroutine corlsd(rs, zta, ec, vcup, vcdn, ecrs, eczta, alfc)
  ! .....-----------------------------------------------------------------
  ! uniform-gas correlation of perdew and wang 1991
  ! .....-----------------------------------------------------------------
  ! input: seitz radius (rs), relative spin polarization (zta)
  ! output: correlation energy per electron (ec),
  ! up- and down-spin potentials (vcup,vcdn),
  ! derivatives of ec wrt rs (ecrs) &zta (eczta).
  ! output: correlation contribution (alfc) to the spin stiffness
  ! .....-----------------------------------------------------------------
  use :: mod_datatypes, only: dp
   use mod_gcor91
  implicit none
  ! .. Scalar Arguments ..
  real (kind=dp) :: alfc, ec, ecrs, eczta, rs, vcdn, vcup, zta
  ! ..
  ! .. Local Scalars ..
  real (kind=dp) :: alfm, alfrsm, comm, ep, eprs, eu, eurs, f, fz, fzz, gam, &
    thrd, thrd4, z4
  ! ..
  ! .. Save statement ..
  save :: gam, fzz, thrd, thrd4
  ! ..
  ! .. Data statements ..
  ! .....-----------------------------------------------------------------
  data gam, fzz/0.5198421e0_dp, 1.709921e0_dp/
  data thrd, thrd4/0.333333333333e0_dp, 1.333333333333e0_dp/
  ! ..
  ! .....-----------------------------------------------------------------
  f = ((1.e0_dp+zta)**thrd4+(1.e0_dp-zta)**thrd4-2.e0_dp)/gam
  call gcor91(0.0310907e0_dp, 0.21370e0_dp, 7.5957e0_dp, 3.5876e0_dp, &
    1.6382e0_dp, 0.49294e0_dp, 1.00e0_dp, rs, eu, eurs)
  call gcor91(0.01554535e0_dp, 0.20548e0_dp, 14.1189e0_dp, 6.1977e0_dp, &
    3.3662e0_dp, 0.62517e0_dp, 1.00e0_dp, rs, ep, eprs)
  call gcor91(0.0168869e0_dp, 0.11125e0_dp, 10.357e0_dp, 3.6231e0_dp, &
    0.88026e0_dp, 0.49671e0_dp, 1.00e0_dp, rs, alfm, alfrsm)
  ! alfm is minus the spin stiffness alfc
  alfc = -alfm
  z4 = zta**4
  ec = eu*(1.e0_dp-f*z4) + ep*f*z4 - alfm*f*(1.e0_dp-z4)/fzz
  ! energy done. now the potential:
  ecrs = eurs*(1.e0_dp-f*z4) + eprs*f*z4 - alfrsm*f*(1.e0_dp-z4)/fzz
  fz = thrd4*((1.e0_dp+zta)**thrd-(1.e0_dp-zta)**thrd)/gam
  eczta = 4.e0_dp*(zta**3)*f*(ep-eu+alfm/fzz) + fz*(z4*ep-z4*eu-(1.e0_dp-z4)* &
    alfm/fzz)
  comm = ec - rs*ecrs/3.e0_dp - zta*eczta
  vcup = comm + eczta
  vcdn = comm - eczta

  return
end subroutine corlsd

end module mod_corlsd
