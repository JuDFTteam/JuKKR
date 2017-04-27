  subroutine check_global_options(itc)
! Generates text describing the global options set by inpsusc.dat
! Checks consistency
  use global

  implicit none

  integer(kind=i4b), intent(in) :: itc
  real(kind=r8b) :: ulen

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
  if (itc == 1) write(*,'("************************************************************")')
! Relativity ahoy
  if (lsoc) then
    if (itc == 1) write(*,'("  SOC   correction  ==>  set atomic options")')
  else
    if (itc == 1) write(*,'("  No SOC correction")')
  end if
! ----------------------------------------------------------------------
! Exchange interactions
  if (ljij) then
    if (itc == 1) write(*,'("  Magnetic couplings Jij  ==>  set atomic options")')
    if (ljijtensor) then
      if (itc == 1) write(*,'("  - Tensor Jij using Ebert formula")')
    else
      if (itc == 1) write(*,'("  - Scalar Jij using Lichtenstein formula")')
    end if
  else
    if (itc == 1) write(*,'("  No magnetic couplings")')
  end if
! ----------------------------------------------------------------------
! External field
  if (lbfield) then
    if (itc == 1) write(*,'("  External B field  ==>  set atomic options")')
  else
    if (itc == 1) write(*,'("  No external B field")')
  end if
! ----------------------------------------------------------------------
! Self-energy model
  if (lsemodel) then
    if (itc == 1) write(*,'("  Model self-energy read in  ==>  provide semodel.dat")')
  else
    if (itc == 1) write(*,'("  No self-energy read in")')
  end if
! ----------------------------------------------------------------------
  if (itc == 1) write(*,'("************************************************************")')
! Susceptibility calculation
  if (lsusc) then
    if (itc == 1) write(*,'("  Susceptibility    ==>  set atomic options")')
    if (itc == 1) write(*,'("  - lmax for susceptibility is",i4)') nlmax0
!   Static or dynamic
    if (ldynsusc) then
      if (itc == 1) write(*,'("  - Dynamic calculation:")')
      if (itc == 1) write(*,'("  - nomega, omegamin, omegamax, domega=",i6,5f10.6)') nomega, omegamin, omegamax, domega
    else
      if (itc == 1) write(*,'("  - Static calculation")')
    end if
!   Kohn-Sham only or enhanced
    if (lenhanced) then
      if (itc == 1) write(*,'("  - Enhanced susceptibility")')
    else
      if (itc == 1) write(*,'("  - Kohn-Sham susceptibility only")')
    end if
!   Which integrals
    if (ldynsusc) then
      if (itc == 1) write(*,'("  - Analytic integral=",l1,", non-analytic integral=",l1)') lanalytic, lnonanalytic
    end if
!   Output
    if (lcartesian) then
      if (itc == 1) write(*,'("  - Susceptibility in cartesian components")')
    else
      if (itc == 1) write(*,'("  - Susceptibility in spin components")')
    end if
!   Shake your Hartree
    if (lkha) then
      if (itc == 1) write(*,'("  Hartree kernel    ==>  set atomic options")')
    else
      if (itc == 1) write(*,'("  No Hartree kernel computed")')
    end if
!   To xc or not to xc, that is the question
    if (lkxc) then
      if (itc == 1) write(*,'("  xc kernel         ==>  set atomic options")')
    else
      if (itc == 1) write(*,'("  No xc kernel computed")')
    end if
!   Sum the spins
    if (lsumrule) then
      if (itc == 1) write(*,'("  Spin sum rule used")')
    else
      if (itc == 1) write(*,'("  Spin sum rule not used ==> DANGER !!!")')
    end if
  else
    if (itc == 1) write(*,'("  No susceptibility calculation")')
!   Set overrides
    lkha = .false.; lkxc = .false.
  end if
! ----------------------------------------------------------------------
  if (itc == 1) write(*,'("************************************************************")')
! Whether to fit the GF
  if (lfit) then
    if (ifit == 1) then
      if (itc == 1) write(*,'("  GF interpolated with n0=",i4)') numd
    else if (ifit == 2) then
      if (itc == 1) write(*,'("  GF fit to rational function with numd, dend=",2i4)') numd, dend
      if (itc == 1) write(*,'("  - Energy argument in rational function shifted by",2f12.6)') eshift
    else
      stop 'Unknown fit/interpolation option'
    end if
!   shift of Re GF
    if (lregf) then
      if (itc == 1) write(*,'("  - Shift of Re GF")')
      if (itc == 1) write(*,'("  - fudge=",f10.6)') fudge
    else
      if (itc == 1) write(*,'("  - No shift of Re GF")')
    end if
  else
    if (itc == 1) write(*,'("  No fitting of GF")')
  end if
  if (itc == 1) write(*,'("************************************************************",/)')
! ----------------------------------------------------------------------
! All done!
  end subroutine check_global_options
