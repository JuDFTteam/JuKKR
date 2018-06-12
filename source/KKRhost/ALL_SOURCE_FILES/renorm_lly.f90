! LLY Lloyd  &
subroutine renorm_lly(cdos_lly, ielast, nspin, natyp, cden, lmaxp1, conc, &
  iestart, ieend, wez, ircut, ipan, ez, zat, rho2ns, r2nef, denef, denefat, &
  espv)
  use :: mod_datatypes, only: dp
  ! Renormalize the valence charge according to Lloyd's formula.
  ! Find renormalization constant per energy, then renormalize
  ! charge/atom/energy, then integrate over energies to find
  ! the renormalized charge/atom. Use it to renormalize the density.
  ! Phivos Mavropoulos, July 2014
  implicit none
  include 'inc.p'
  integer :: lmaxd1
  parameter (lmaxd1=lmaxd+1)
  integer :: lmpotd
  parameter (lmpotd=(lpotd+1)**2)
  integer :: npotd
  parameter (npotd=(2*(krel+korbit)+(1-(krel+korbit))*nspind)*natypd)
  ! Concentration (for cpa)
  integer :: lmaxp1, natyp, nspin
  integer :: iestart, ieend, ielast ! Non-renormalized density per atom
                                    ! (density=-cden/pi)
  integer :: ircut(0:ipand, natypd), ipan(natypd) ! DOS according to Lloyd's
                                                  ! formula
  real (kind=dp) :: conc(natypd)   ! Input/Output:
  complex (kind=dp) :: cden(0:lmaxd1, ielast, npotd) ! Internal:
  complex (kind=dp) :: cdos_lly(iemxd, nspind) ! 1: charge renormalization per
                                               ! atom (energy-integrated)
  complex (kind=dp) :: wez(iemxd), ez(iemxd)
  real (kind=dp) :: zat(natypd)
  ! 2: same for spin moment
  real (kind=dp) :: rho2ns(irmd, lmpotd, natypd, 2)
  real (kind=dp) :: r2nef(irmd, lmpotd, natypd, 2)
  real (kind=dp) :: denef, denefat(natypd)
  real (kind=dp) :: espv(0:lmaxd1, npotd)
  ! Density from local summation
  integer :: ll, ie, i1, ispin, ipot, spindegen, irc1, signsp, idim
  real (kind=dp) :: renorm_at(natypd, 2) ! and from Lloyd's formula
  ! Renormalization constant for charge and spin density
  complex (kind=dp) :: cdos_loc(iemxd, (1+krel)*nspind) ! Atomic charge per
                                                        ! spin (local
                                                        ! summation and
                                                        ! renormalized)
  complex (kind=dp) :: cdos_locvc(iemxd, (1+krel)*nspind)
  real (kind=dp) :: cren(iemxd, 2)
  real (kind=dp) :: charge(natypd, 2), charge_lly(natypd, 2)
  complex (kind=dp) :: chadd(iemxd, natypd, nspind), cdos_add ! Integration
                                                              ! step for
                                                              ! charge/atom/spin
  complex (kind=dp) :: qlly(2), qstar(2)
  real (kind=dp) :: sum0(2), sum1(2)
  complex (kind=dp) :: czero
  real (kind=dp) :: pi
  logical :: opt, test
  external :: opt, test
  ! Spin degeneracy, 2 if nspin=1, 1 if nspin=2

  czero = (0.e0_dp, 0.e0_dp)
  pi = 4.e0_dp*atan(1.e0_dp)

  spindegen = 3 - nspin            ! First find renormalization factor per
                                   ! energy and atomic charges
  ! Factor 1/pi included in Wez
  cren(:, :) = 0e0_dp
  renorm_at(:, :) = 1.e0_dp
  charge_lly(:, :) = 0.e0_dp
  charge(:, :) = 0.e0_dp
  qlly(:) = czero
  qstar(:) = czero
  ! Complex charge
  ! I1=1,NATYP
  cdos_loc = czero
  cdos_locvc = czero
  chadd = czero
  do ie = iestart, ieend
    do ispin = 1, nspin
      do i1 = 1, natyp
        ipot = (i1-1)*nspin + ispin
        cdos_add = czero
        do ll = 0, lmaxp1
          cdos_add = cdos_add + conc(i1)*cden(ll, ie, ipot) ! ISPIN = 1,NSPIN
        end do
        if (zat(i1)>1e-06_dp) then
          cdos_loc(ie, ispin) = cdos_loc(ie, ispin) + cdos_add
        else
          cdos_locvc(ie, ispin) = cdos_locvc(ie, ispin) + cdos_add
        end if
        chadd(ie, i1, ispin) = wez(ie)*cdos_add ! IE = IESTART,IEEND
        charge(i1, ispin) = charge(i1, ispin) + aimag(chadd(ie,i1,ispin))/real &
          (nspin, kind=dp)
      end do                       ! Now the locally-summed charge/energy is
                                   ! in cdos_loc, charge/energy/atom in chadd
      cdos_loc(ie, ispin) = -cdos_loc(ie, ispin)/pi
      cdos_locvc(ie, ispin) = -cdos_locvc(ie, ispin)/pi
    end do                         ! Renormalization factor per energy:
  end do                           ! Apply to DOS of each atom:
  ! ISPIN = 1,NSPIN
  if (.not. opt('NEWSOSOL')) then
    do ie = iestart, ieend
      do ispin = 1, nspin
        ! IE = IESTART,IEEND
        cren(ie, ispin) = aimag((cdos_lly(ie,ispin)-cdos_locvc(ie, &
          ispin))*wez(ie))/aimag(cdos_loc(ie,ispin)*wez(ie))
        ! Renormalization factor per energy:
        do i1 = 1, natypd
          if (zat(i1)>1e-06_dp) then
            charge_lly(i1, ispin) = charge_lly(i1, ispin) + &
              cren(ie, ispin)*aimag(chadd(ie,i1,ispin))/real(nspin, kind=dp)
          else
            charge_lly(i1, ispin) = charge_lly(i1, ispin) + &
              aimag(chadd(ie,i1,ispin))/real(nspin, kind=dp)
          end if
        end do
      end do                       ! Apply to DOS of each atom:
    end do
  else
    do ie = iestart, ieend
      ! IE = IESTART,IEEND
      cren(ie, 1) = aimag((cdos_lly(ie,1)-cdos_locvc(ie,1)-cdos_locvc(ie, &
        2))*wez(ie))/aimag((cdos_loc(ie,1)+cdos_loc(ie,2))*wez(ie))
      ! add term from sum from l>lmax to infinity
      do ispin = 1, nspin
        do i1 = 1, natypd
          if (zat(i1)>1e-06_dp) then
            charge_lly(i1, ispin) = charge_lly(i1, ispin) + &
              cren(ie, 1)*aimag(chadd(ie,i1,ispin))/real(nspin, kind=dp)
          else
            charge_lly(i1, ispin) = charge_lly(i1, ispin) + &
              aimag(chadd(ie,i1,ispin))/real(nspin, kind=dp)
          end if
        end do
      end do
    end do                         ! DO I1=1,NATYPD
  end if
  ! DO ISPIN=1,NSPIN
  ! CHARGE_LLY(I1,ISPIN)=CHARGE_LLY(I1,ISPIN)-DIMAG(CDOS2(I1))
  ! ENDDO
  ! ENDDO


  ! Now apply renormalization to energy-integrated density
  ! If spins are coupled, then only charge density
  if (nspin==1 .or. opt('NEWSOSOL')) cren(:, 2) = cren(:, 1)
  ! Index of outmost radial point


  if (nspin==1) then
    do i1 = 1, natyp
      if (charge(i1,1)>0) then
        renorm_at(i1, 1) = charge_lly(i1, 1)/charge(i1, 1)
      else
        renorm_at(i1, 1) = 1.0e0_dp
      end if
      renorm_at(i1, 2) = renorm_at(i1, 1)
      irc1 = ircut(ipan(i1), i1)
      rho2ns(1:irc1, 1:lmpotd, i1, 1) = rho2ns(1:irc1, 1:lmpotd, i1, 1)* &
        renorm_at(i1, 1)
      r2nef(1:irc1, 1:lmpotd, i1, 1) = r2nef(1:irc1, 1:lmpotd, i1, 1)* &
        renorm_at(i1, 1)
    end do
    ! First decouple charge and spin density to the density of each channels
  else
    ! Index of outmost radial point
    ! Second merge density of each channels to charge and spin density
    idim = irmd*lmpotd*natypd
    call daxpy(idim, 1.0e0_dp, rho2ns(1,1,1,1), 1, rho2ns(1,1,1,2), 1)
    call dscal(idim, 0.5e0_dp, rho2ns(1,1,1,2), 1)
    call daxpy(idim, -1.0e0_dp, rho2ns(1,1,1,2), 1, rho2ns(1,1,1,1), 1)

    call daxpy(idim, 1.0e0_dp, r2nef(1,1,1,1), 1, r2nef(1,1,1,2), 1)
    call dscal(idim, 0.5e0_dp, r2nef(1,1,1,2), 1)
    call daxpy(idim, -1.0e0_dp, r2nef(1,1,1,2), 1, r2nef(1,1,1,1), 1)
    do i1 = 1, natyp
      irc1 = ircut(ipan(i1), i1)
      do ispin = 1, nspin
        renorm_at(i1, ispin) = charge_lly(i1, ispin)/charge(i1, ispin)
        rho2ns(1:irc1, 1:lmpotd, i1, ispin) = rho2ns(1:irc1, 1:lmpotd, i1, &
          ispin)*renorm_at(i1, ispin)
        r2nef(1:irc1, 1:lmpotd, i1, ispin) = r2nef(1:irc1, 1:lmpotd, i1, &
          ispin)*renorm_at(i1, ispin)
      end do
    end do
    ! calculate density at Fermi level
    call dscal(idim, 2.0e0_dp, rho2ns(1,1,1,1), 1)
    call daxpy(idim, -0.5e0_dp, rho2ns(1,1,1,1), 1, rho2ns(1,1,1,2), 1)
    call daxpy(idim, 1.0e0_dp, rho2ns(1,1,1,2), 1, rho2ns(1,1,1,1), 1)
    ! LL
    call dscal(idim, 2.0e0_dp, r2nef(1,1,1,1), 1)
    call daxpy(idim, -0.5e0_dp, r2nef(1,1,1,1), 1, r2nef(1,1,1,2), 1)
    call daxpy(idim, 1.0e0_dp, r2nef(1,1,1,2), 1, r2nef(1,1,1,1), 1)
  end if
  ! ISPIN
  ! I1
  denef = 0e0_dp
  do i1 = 1, natyp
    denefat(i1) = 0e0_dp
    do ispin = 1, nspin
      ipot = (i1-1)*nspin + ispin
      do ll = 0, lmaxp1
        if (zat(i1)>1e-06_dp) then
          denefat(i1) = denefat(i1) - 2.0e0_dp*conc(i1)*cren(ielast, ispin)* &
            aimag(cden(ll,ielast,ipot))/pi/real(nspin, kind=dp)
        else
          denefat(i1) = denefat(i1) - 2.0e0_dp*conc(i1)*aimag(cden(ll,ielast, &
            ipot))/pi/real(nspin, kind=dp)
        end if
        espv(ll, ipot) = 0e0_dp
        if (zat(i1)>1e-06_dp) then
          do ie = 1, ielast
            espv(ll, ipot) = espv(ll, ipot) + cren(ie, ispin)*aimag(ez(ie)* &
              cden(ll,ie,ipot)*wez(ie)/real(nspin,kind=dp))
          end do
        else
          do ie = 1, ielast
            espv(ll, ipot) = espv(ll, ipot) + aimag(ez(ie)*cden(ll,ie,ipot)* &
              wez(ie)/real(nspin,kind=dp))
          end do
        end if
      end do                       ! Write out renormalization factors
    end do
    denef = denef + denefat(i1)
  end do
  ! -1,+1 for spin down,up (ispin=1,2)
  write (1337, *) 'Information on renormalization by Lloyds formula'
  write (1337, *) 'RENORM_LLY: Complex renormalization factor per energy:'
  write (1337, fmt='(A5,2A32)') 'IE', 'Spin 1 (down)           ', &
    'Spin 2 (up)           '
  do ie = iestart, ieend
    write (1337, fmt='(I5,4F16.12)') ie, (cren(ie,ispin), ispin=1, nspin)
  end do
  write (1337, *) 'RENORM_LLY: renormalization factor per atom:'
  write (1337, fmt='(A5,2A16)') 'IAT', 'Spin down', 'Spin up'
  do i1 = 1, natyp
    write (1337, fmt='(I5,2E17.9)') i1, (renorm_at(i1,ispin), ispin=1, nspin)
  end do
  write (1337, *) 'RENORM_LLY: Renormalized charge per atom:'
  write (1337, fmt='(A5,2A16)') 'IAT', 'Spin down', 'Spin up'
  do i1 = 1, natyp
    write (1337, fmt='(I5,2F16.12)') i1, (charge_lly(i1,ispin), ispin=1, nspin &
      )
  end do
  sum0(:) = 0.e0_dp
  sum1(:) = 0.e0_dp
  do ispin = 1, nspin
    signsp = 2*ispin - 3
    if (nspin==1) signsp = 1
    do i1 = 1, natyp
      sum0(ispin) = sum0(ispin) + signsp*conc(i1)*charge(i1, ispin)
      sum1(ispin) = sum1(ispin) + signsp*conc(i1)*charge_lly(i1, ispin)

    end do
  end do
  write (1337, fmt='(A45,2E17.9)') &
    'RENORM_LLY: Locally summed charge and moment:', &
    (sum0(ispin), ispin=1, nspin)
  write (1337, fmt='(A45,2E17.9)') &
    'RENORM_LLY: Renormalized charge and moment:  ', &
    (sum1(ispin), ispin=1, nspin)
  write (1337, fmt='(A50,2E17.9)') &
    'RENORM_LLY: Renormalization factor of total charge:', sum1(1)/sum0(1)
  ! LLY Lloyd  &
  ! Renormalize the valence charge according to Lloyd's formula.
  ! Find renormalization constant per energy, then renormalize
end subroutine renorm_lly
