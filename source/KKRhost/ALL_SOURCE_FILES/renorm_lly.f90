subroutine renorm_lly( & ! LLY Lloyd  &
  cdos_lly, ielast, nspin, natyp, cden, lmaxp1, conc, iestart, ieend, wez, &
  ircut, ipan, ez, zat, rho2ns, r2nef, denef, denefat, espv)
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
  integer :: lmaxp1, natyp, nspin ! Non-renormalized density per atom (density=-cden/pi)
  integer :: iestart, ieend, ielast ! DOS according to Lloyd's formula
  integer :: ircut(0:ipand, natypd), ipan(natypd) ! Input/Output:
  double precision :: conc(natypd) ! Internal:
  double complex :: cden(0:lmaxd1, ielast, npotd) ! 1: charge renormalization per atom (energy-integrated)
  double complex :: cdos_lly(iemxd, nspind) ! 2: same for spin moment
  double complex :: wez(iemxd), ez(iemxd)
  double precision :: zat(natypd)
!  Density from local summation
  double precision :: rho2ns(irmd, lmpotd, natypd, 2)
  double precision :: r2nef(irmd, lmpotd, natypd, 2)
  double precision :: denef, denefat(natypd)
  double precision :: espv(0:lmaxd1, npotd)
!      and from Lloyd's formula
  integer :: ll, ie, i1, ispin, ipot, spindegen, irc1, signsp, idim
  double precision :: renorm_at(natypd, 2) ! Renormalization constant for charge and spin density
! Atomic charge per spin (local summation and renormalized)
  double complex :: cdos_loc(iemxd, (1+krel)*nspind) ! Integration step for charge/atom/spin
  double complex :: cdos_locvc(iemxd, (1+krel)*nspind) 
  double precision :: cren(iemxd, 2) 
  double precision :: charge(natypd, 2), charge_lly(natypd, 2) 
  double complex :: chadd(iemxd, natypd, nspind), cdos_add ! Spin degeneracy, 2 if nspin=1, 1 if nspin=2
  double complex :: qlly(2), qstar(2)
  double precision :: sum0(2), sum1(2)
  double complex :: czero
  double precision :: pi
  logical :: opt, test
  external :: opt, test


  czero = (0.d0, 0.d0)
  pi = 4.d0*datan(1.d0)
! First find renormalization factor per energy and atomic charges
  spindegen = 3 - nspin ! Factor 1/pi included in Wez
! Complex charge
  cren(:, :) = 0d0
  renorm_at(:, :) = 1.d0
  charge_lly(:, :) = 0.d0
  charge(:, :) = 0.d0
  qlly(:) = czero
  qstar(:) = czero
! I1=1,NATYP
! ISPIN = 1,NSPIN
  cdos_loc = czero
  cdos_locvc = czero
  chadd = czero
  do ie = iestart, ieend
    do ispin = 1, nspin
      do i1 = 1, natyp
        ipot = (i1-1)*nspin + ispin
        cdos_add = czero
        do ll = 0, lmaxp1
          cdos_add = cdos_add + conc(i1)*cden(ll, ie, ipot) ! IE = IESTART,IEEND
        end do
        if (zat(i1)>1d-06) then
          cdos_loc(ie, ispin) = cdos_loc(ie, ispin) + cdos_add
        else
          cdos_locvc(ie, ispin) = cdos_locvc(ie, ispin) + cdos_add
        end if
        chadd(ie, i1, ispin) = wez(ie)*cdos_add ! Now the locally-summed charge/energy is in cdos_loc, charge/energy/atom in chadd
        charge(i1, ispin) = charge(i1, ispin) + dimag(chadd(ie,i1,ispin))/dble &
          (nspin)
      end do ! Renormalization factor per energy:
      cdos_loc(ie, ispin) = -cdos_loc(ie, ispin)/pi
      cdos_locvc(ie, ispin) = -cdos_locvc(ie, ispin)/pi
    end do ! Apply to DOS of each atom:
  end do ! ISPIN = 1,NSPIN
! IE = IESTART,IEEND
  if (.not. opt('NEWSOSOL')) then
    do ie = iestart, ieend
      do ispin = 1, nspin
! Renormalization factor per energy:
        cren(ie, ispin) = dimag((cdos_lly(ie,ispin)-cdos_locvc(ie, &
          ispin))*wez(ie))/dimag(cdos_loc(ie,ispin)*wez(ie))
! Apply to DOS of each atom:
        do i1 = 1, natypd
          if (zat(i1)>1d-06) then
            charge_lly(i1, ispin) = charge_lly(i1, ispin) + &
              cren(ie, ispin)*dimag(chadd(ie,i1,ispin))/dble(nspin)
          else
            charge_lly(i1, ispin) = charge_lly(i1, ispin) + &
              dimag(chadd(ie,i1,ispin))/dble(nspin)
          end if
        end do
      end do ! IE = IESTART,IEEND
    end do 
  else
    do ie = iestart, ieend
! add term from sum from l>lmax to infinity
      cren(ie, 1) = dimag((cdos_lly(ie,1)-cdos_locvc(ie,1)-cdos_locvc(ie, &
        2))*wez(ie))/dimag((cdos_loc(ie,1)+cdos_loc(ie,2))*wez(ie))
!            DO I1=1,NATYPD
      do ispin = 1, nspin
        do i1 = 1, natypd
          if (zat(i1)>1d-06) then
            charge_lly(i1, ispin) = charge_lly(i1, ispin) + &
              cren(ie, 1)*dimag(chadd(ie,i1,ispin))/dble(nspin)
          else
            charge_lly(i1, ispin) = charge_lly(i1, ispin) + &
              dimag(chadd(ie,i1,ispin))/dble(nspin)
          end if
        end do
      end do
    end do !             DO ISPIN=1,NSPIN
  end if
!             CHARGE_LLY(I1,ISPIN)=CHARGE_LLY(I1,ISPIN)-DIMAG(CDOS2(I1))
!             ENDDO
!            ENDDO


! Now apply renormalization to energy-integrated density
! If spins are coupled, then only charge density
! Index of outmost radial point
  if (nspin==1 .or. opt('NEWSOSOL')) cren(:, 2) = cren(:, 1)


! First decouple charge and spin density to the density of each channels
  if (nspin==1) then
    do i1 = 1, natyp
      if (charge(i1,1)>0) then
        renorm_at(i1, 1) = charge_lly(i1, 1)/charge(i1, 1)
      else
        renorm_at(i1, 1) = 1.0d0
      end if
      renorm_at(i1, 2) = renorm_at(i1, 1)
      irc1 = ircut(ipan(i1), i1) 
      rho2ns(1:irc1, 1:lmpotd, i1, 1) = rho2ns(1:irc1, 1:lmpotd, i1, 1)* &
        renorm_at(i1, 1)
      r2nef(1:irc1, 1:lmpotd, i1, 1) = r2nef(1:irc1, 1:lmpotd, i1, 1)* &
        renorm_at(i1, 1)
    end do
! Index of outmost radial point
  else
! Second merge density of each channels to charge and spin density

    idim = irmd*lmpotd*natypd
    call daxpy(idim, 1.0d0, rho2ns(1,1,1,1), 1, rho2ns(1,1,1,2), 1)
    call dscal(idim, 0.5d0, rho2ns(1,1,1,2), 1)
    call daxpy(idim, -1.0d0, rho2ns(1,1,1,2), 1, rho2ns(1,1,1,1), 1)

    call daxpy(idim, 1.0d0, r2nef(1,1,1,1), 1, r2nef(1,1,1,2), 1)
    call dscal(idim, 0.5d0, r2nef(1,1,1,2), 1)
    call daxpy(idim, -1.0d0, r2nef(1,1,1,2), 1, r2nef(1,1,1,1), 1)
    do i1 = 1, natyp
      irc1 = ircut(ipan(i1), i1) ! calculate density at Fermi level
      do ispin = 1, nspin
        renorm_at(i1, ispin) = charge_lly(i1, ispin)/charge(i1, ispin)
        rho2ns(1:irc1, 1:lmpotd, i1, ispin) = rho2ns(1:irc1, 1:lmpotd, i1, &
          ispin)*renorm_at(i1, ispin)
        r2nef(1:irc1, 1:lmpotd, i1, ispin) = r2nef(1:irc1, 1:lmpotd, i1, &
          ispin)*renorm_at(i1, ispin)
      end do
    end do
! LL
    call dscal(idim, 2.0d0, rho2ns(1,1,1,1), 1)
    call daxpy(idim, -0.5d0, rho2ns(1,1,1,1), 1, rho2ns(1,1,1,2), 1)
    call daxpy(idim, 1.0d0, rho2ns(1,1,1,2), 1, rho2ns(1,1,1,1), 1)
! ISPIN
    call dscal(idim, 2.0d0, r2nef(1,1,1,1), 1)
    call daxpy(idim, -0.5d0, r2nef(1,1,1,1), 1, r2nef(1,1,1,2), 1)
    call daxpy(idim, 1.0d0, r2nef(1,1,1,2), 1, r2nef(1,1,1,1), 1)
  end if
! I1
! Write out renormalization factors
  denef = 0d0
  do i1 = 1, natyp
    denefat(i1) = 0d0
    do ispin = 1, nspin
      ipot = (i1-1)*nspin + ispin
      do ll = 0, lmaxp1
        if (zat(i1)>1d-06) then
          denefat(i1) = denefat(i1) - 2.0d0*conc(i1)*cren(ielast, ispin)*dimag &
            (cden(ll,ielast,ipot))/pi/dble(nspin)
        else
          denefat(i1) = denefat(i1) - 2.0d0*conc(i1)*dimag(cden(ll,ielast,ipot &
            ))/pi/dble(nspin)
        end if
        espv(ll, ipot) = 0d0
        if (zat(i1)>1d-06) then
          do ie = 1, ielast
            espv(ll, ipot) = espv(ll, ipot) + cren(ie, ispin)*dimag(ez(ie)* &
              cden(ll,ie,ipot)*wez(ie)/dble(nspin))
          end do
        else
          do ie = 1, ielast
            espv(ll, ipot) = espv(ll, ipot) + dimag(ez(ie)*cden(ll,ie,ipot)* &
              wez(ie)/dble(nspin))
          end do
        end if
      end do ! -1,+1 for spin down,up (ispin=1,2)
    end do 
    denef = denef + denefat(i1)
  end do 

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
  sum0(:) = 0.d0
  sum1(:) = 0.d0
  do ispin = 1, nspin
    signsp = 2*ispin - 3 
    if (nspin==1) signsp = 1
    do i1 = 1, natyp
      sum0(ispin) = sum0(ispin) + signsp*conc(i1)*charge(i1, ispin)
      sum1(ispin) = sum1(ispin) + signsp*conc(i1)*charge_lly(i1, ispin)
! LLY Lloyd  &
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
! Renormalize the valence charge according to Lloyd's formula.
! Find renormalization constant per energy, then renormalize
! charge/atom/energy, then integrate over energies to find
end subroutine
