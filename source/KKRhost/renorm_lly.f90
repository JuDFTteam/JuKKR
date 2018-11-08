!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Renormalize the valence charge according to Lloyd's formula.
!> Author: Phivos Mavropoulos
!> Renormalize the valence charge according to Lloyd's formula. 
!> Find renormalization constant per energy, then renormalize charge/atom/energy, 
!> then integrate over energies to find the renormalized charge/atom. 
!> Use it to renormalize the density.
!------------------------------------------------------------------------------------
module mod_renorm_lly

contains

  !-------------------------------------------------------------------------------
  !> Summary: Renormalize the valence charge according to Lloyd's formula.
  !> Author: Phivos Mavropoulos
  !> Category: physical-observables, KKRhost
  !> Deprecated: False 
  !> Renormalize the valence charge according to Lloyd's formula. 
  !> Find renormalization constant per energy, then renormalize charge/atom/energy, 
  !> then integrate over energies to find the renormalized charge/atom. 
  !> Use it to renormalize the density.
  !-------------------------------------------------------------------------------
  subroutine renorm_lly(cdos_lly, ielast, nspin, natyp, cden, lmaxp1, conc, iestart, ieend, wez, ircut, ipan, ez, zat, rho2ns, r2nef, denef, denefat, espv)

    use :: mod_datatypes, only: dp
    use :: mod_runoptions, only: use_Chebychev_solver, set_cheby_nosoc
    use :: mod_constants, only: czero, pi
    use :: global_variables, only: ipand, natypd, lmaxd, npotd, iemxd, irmd, lmpotd, krel, nspind 
    implicit none
   
    integer :: lmaxp1, natyp, nspin
    integer :: iestart, ieend, ielast

    integer :: ircut(0:ipand, natypd), ipan(natypd)

    real (kind=dp) :: conc(natypd) !! Concentration (for cpa)
    complex (kind=dp) :: cden(0:(lmaxd+1), ielast, npotd) !! Non-renormalized density per atom (density=-cden/pi)
    complex (kind=dp) :: cdos_lly(ielast, nspin) !! DOS according to Lloyd's formula
    complex (kind=dp) :: wez(iemxd), ez(iemxd)
    real (kind=dp) :: zat(natypd)
! Input/Output:
    real (kind=dp) :: rho2ns(irmd, lmpotd, natypd, 2)
    real (kind=dp) :: r2nef(irmd, lmpotd, natypd, 2)
    real (kind=dp) :: denef, denefat(natypd)
    real (kind=dp) :: espv(0:(lmaxd+1), npotd)
! Internal:
    integer :: ll, ie, i1, ispin, ipot, spindegen, irc1, signsp, idim
    real (kind=dp) :: renorm_at(natypd, 2) !! 1: charge renormalization per atom (energy-integrated); 2: same for spin moment
    complex (kind=dp) :: cdos_loc(ielast, (1+krel)*nspind) !! Density from local summation
    complex (kind=dp) :: cdos_locvc(ielast, (1+krel)*nspind) !! Density from Lloyd's formula
    real (kind=dp) :: cren(ielast, 2) !! Renormalization constant for charge and spin density
    real (kind=dp) :: charge(natypd, 2), charge_lly(natypd, 2) !! Atomic charge per spin (local summation and renormalized)
    complex (kind=dp) :: chadd(ielast, natypd, nspind), cdos_add !! Integration step for charge/atom/spin
    complex (kind=dp) :: qlly(2), qstar(2)
    real (kind=dp) :: sum0(2), sum1(2)
    logical, external :: opt, test


    ! Spin degeneracy, 2 if nspin=1, 1 if nspin=2
    spindegen = 3 - nspin

    cren(:, :) = 0e0_dp
    renorm_at(:, :) = 1.e0_dp
    charge_lly(:, :) = 0.e0_dp
    charge(:, :) = 0.e0_dp
    qlly(:) = czero
    qstar(:) = czero

    ! First find renormalization factor per energy and atomic charges
    cdos_loc = czero
    cdos_locvc = czero
    chadd = czero
    do ie = iestart, ieend
      do ispin = 1, nspin
        do i1 = 1, natyp
          ipot = (i1-1)*nspin + ispin
          cdos_add = czero
          do ll = 0, lmaxp1
            ! Factor 1/pi included in Wez
            cdos_add = cdos_add + conc(i1)*cden(ll, ie, ipot)
          end do
          if (zat(i1)>1e-06_dp) then
            cdos_loc(ie, ispin) = cdos_loc(ie, ispin) + cdos_add
          else
            cdos_locvc(ie, ispin) = cdos_locvc(ie, ispin) + cdos_add
          end if
          ! Complex charge
          chadd(ie, i1, ispin) = wez(ie)*cdos_add
          charge(i1, ispin) = charge(i1, ispin) + aimag(chadd(ie,i1,ispin))/real(nspin, kind=dp)
        end do
        cdos_loc(ie, ispin) = -cdos_loc(ie, ispin)/pi
        cdos_locvc(ie, ispin) = -cdos_locvc(ie, ispin)/pi
      end do
    end do
    ! Now the locally-summed charge/energy is in cdos_loc, charge/energy/atom in chadd
    if (.not. use_Chebychev_solver .or. set_cheby_nosoc) then
      do ie = iestart, ieend
        do ispin = 1, nspin
          ! Renormalization factor per energy:
          cren(ie, ispin) = aimag((cdos_lly(ie,ispin)-cdos_locvc(ie,ispin))*wez(ie))/aimag(cdos_loc(ie,ispin)*wez(ie))
          ! Apply to DOS of each atom:
          do i1 = 1, natypd
            if (zat(i1)>1e-06_dp) then
              charge_lly(i1, ispin) = charge_lly(i1, ispin) + cren(ie, ispin)*aimag(chadd(ie,i1,ispin))/real(nspin, kind=dp)
            else
              charge_lly(i1, ispin) = charge_lly(i1, ispin) + aimag(chadd(ie,i1,ispin))/real(nspin, kind=dp)
            end if
          end do
        end do
      end do
    else
      do ie = iestart, ieend
        ! Renormalization factor per energy:
        cren(ie, 1) = aimag((cdos_lly(ie,1)-cdos_locvc(ie,1)-cdos_locvc(ie,2))*wez(ie))/aimag((cdos_loc(ie,1)+cdos_loc(ie,2))*wez(ie))
        ! Apply to DOS of each atom:
        do ispin = 1, nspin
          do i1 = 1, natypd
            if (zat(i1)>1e-06_dp) then
              charge_lly(i1, ispin) = charge_lly(i1, ispin) + cren(ie, 1)*aimag(chadd(ie,i1,ispin))/real(nspin, kind=dp)
            else
              charge_lly(i1, ispin) = charge_lly(i1, ispin) + aimag(chadd(ie,i1,ispin))/real(nspin, kind=dp)
            end if
          end do
        end do
      end do
    end if

    ! add term from sum from l>lmax to infinity
    !   DO I1=1,NATYPD
    !    DO ISPIN=1,NSPIN
    !    CHARGE_LLY(I1,ISPIN)=CHARGE_LLY(I1,ISPIN)-DIMAG(CDOS2(I1))
    !    ENDDO
    !   ENDDO

    if (nspin==1 .or. (use_Chebychev_solver .and. .not. set_cheby_nosoc) ) cren(:, 2) = cren(:, 1)


    ! Now apply renormalization to energy-integrated density
    ! If spins are coupled, then only charge density
    if (nspin==1) then
      do i1 = 1, natyp
        if (charge(i1,1)>0) then
          renorm_at(i1, 1) = charge_lly(i1, 1)/charge(i1, 1)
        else
          renorm_at(i1, 1) = 1.0e0_dp
        end if
        renorm_at(i1, 2) = renorm_at(i1, 1)
        irc1 = ircut(ipan(i1), i1) ! Index of outmost radial point
        rho2ns(1:irc1, 1:lmpotd, i1, 1) = rho2ns(1:irc1, 1:lmpotd, i1, 1)*renorm_at(i1, 1)
        r2nef(1:irc1, 1:lmpotd, i1, 1) = r2nef(1:irc1, 1:lmpotd, i1, 1)*renorm_at(i1, 1)
      end do
    else
      ! First decouple charge and spin density to the density of each channels
      idim = irmd*lmpotd*natypd
      call daxpy(idim, 1.0e0_dp, rho2ns(1,1,1,1), 1, rho2ns(1,1,1,2), 1)
      call dscal(idim, 0.5e0_dp, rho2ns(1,1,1,2), 1)
      call daxpy(idim, -1.0e0_dp, rho2ns(1,1,1,2), 1, rho2ns(1,1,1,1), 1)

      call daxpy(idim, 1.0e0_dp, r2nef(1,1,1,1), 1, r2nef(1,1,1,2), 1)
      call dscal(idim, 0.5e0_dp, r2nef(1,1,1,2), 1)
      call daxpy(idim, -1.0e0_dp, r2nef(1,1,1,2), 1, r2nef(1,1,1,1), 1)
      do i1 = 1, natyp
        irc1 = ircut(ipan(i1), i1) ! Index of outmost radial point
        do ispin = 1, nspin
          renorm_at(i1, ispin) = charge_lly(i1, ispin)/charge(i1, ispin)
          rho2ns(1:irc1, 1:lmpotd, i1, ispin) = rho2ns(1:irc1, 1:lmpotd, i1, ispin)*renorm_at(i1, ispin)
          r2nef(1:irc1, 1:lmpotd, i1, ispin) = r2nef(1:irc1, 1:lmpotd, i1, ispin)*renorm_at(i1, ispin)
        end do
      end do
      ! Second merge density of each channels to charge and spin density
      call dscal(idim, 2.0e0_dp, rho2ns(1,1,1,1), 1)
      call daxpy(idim, -0.5e0_dp, rho2ns(1,1,1,1), 1, rho2ns(1,1,1,2), 1)
      call daxpy(idim, 1.0e0_dp, rho2ns(1,1,1,2), 1, rho2ns(1,1,1,1), 1)

      call dscal(idim, 2.0e0_dp, r2nef(1,1,1,1), 1)
      call daxpy(idim, -0.5e0_dp, r2nef(1,1,1,1), 1, r2nef(1,1,1,2), 1)
      call daxpy(idim, 1.0e0_dp, r2nef(1,1,1,2), 1, r2nef(1,1,1,1), 1)
    end if

    ! calculate density at Fermi level
    denef = 0e0_dp
    do i1 = 1, natyp
      denefat(i1) = 0e0_dp
      do ispin = 1, nspin
        ipot = (i1-1)*nspin + ispin
        do ll = 0, lmaxp1
          if (zat(i1)>1e-06_dp) then
            denefat(i1) = denefat(i1) - 2.0e0_dp*conc(i1)*cren(ielast, ispin)*aimag(cden(ll,ielast,ipot))/pi/real(nspin, kind=dp)
          else
            denefat(i1) = denefat(i1) - 2.0e0_dp*conc(i1)*aimag(cden(ll,ielast,ipot))/pi/real(nspin, kind=dp)
          end if
          espv(ll, ipot) = 0e0_dp
          if (zat(i1)>1e-06_dp) then
            do ie = 1, ielast
              espv(ll, ipot) = espv(ll, ipot) + cren(ie, ispin)*aimag(ez(ie)*cden(ll,ie,ipot)*wez(ie)/real(nspin,kind=dp))
            end do
          else
            do ie = 1, ielast
              espv(ll, ipot) = espv(ll, ipot) + aimag(ez(ie)*cden(ll,ie,ipot)*wez(ie)/real(nspin,kind=dp))
            end do
          end if
        end do
      end do
      denef = denef + denefat(i1)
    end do

    ! Write out renormalization factors
    write (1337, *) 'Information on renormalization by Lloyds formula'
    write (1337, *) 'RENORM_LLY: Complex renormalization factor per energy:'
    write (1337, fmt='(A5,2A32)') 'IE', 'Spin 1 (down)           ', 'Spin 2 (up)           '
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
      write (1337, fmt='(I5,2F16.12)') i1, (charge_lly(i1,ispin), ispin=1, nspin)
    end do
    sum0(:) = 0.e0_dp
    sum1(:) = 0.e0_dp
    do ispin = 1, nspin
      signsp = 2*ispin - 3 ! -1,+1 for spin down,up (ispin=1,2)
      if (nspin==1) signsp = 1
      do i1 = 1, natyp
        sum0(ispin) = sum0(ispin) + signsp*conc(i1)*charge(i1, ispin)
        sum1(ispin) = sum1(ispin) + signsp*conc(i1)*charge_lly(i1, ispin)
      end do
    end do
    write (1337, fmt='(A45,2E17.9)') 'RENORM_LLY: Locally summed charge and moment:', (sum0(ispin), ispin=1, nspin)
    write (1337, fmt='(A45,2E17.9)') 'RENORM_LLY: Renormalized charge and moment:  ', (sum1(ispin), ispin=1, nspin)
    write (1337, fmt='(A50,2E17.9)') 'RENORM_LLY: Renormalization factor of total charge:', sum1(1)/sum0(1)

  end subroutine renorm_lly

end module mod_renorm_lly
