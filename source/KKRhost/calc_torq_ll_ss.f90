!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_calc_torq_ll_ss

contains

  !-------------------------------------------------------------------------------
  !> Summary: Torque matrix for PKKprime code
  !> Author: G. Geranton
  !> Date: September 2014
  !> Category: KKRhost, physical-observables
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !> 
  !> This subroutine computes a matrix that is the basis for constructing
  !> the KKR representation of the torque operator. It is adapted from the
  !> CALC_RHO_LL_SS subroutine, but the spin dependent part, i.e., the exhange
  !> field, replaces the shape function in the integration.
  !>
  !> Guillaume Geranton, September 2014
  !-------------------------------------------------------------------------------
  subroutine calc_torq_ll_ss(lmmax0d, rll, ircut, ipan, icell, cleb, icleb, iend, lmsp, irws, drdi, dens, visp, nspin, iatom, vins, irmin)

    use :: mod_runoptions, only: torque_operator_onlyMT, torque_operator_onlySph 
    use :: global_variables, only: irmd, lmpotd, irmind, ipand, natypd, ncleb
    use :: mod_datatypes, only: dp
    use :: mod_csimpk, only: csimpk

    implicit none

    ! .. Array Arguments ..
    ! non-sph. eigen states of single pot  &
    integer :: iend, lmmax0d, irws, nspin, iatom, irmin ! derivative dr/di  &
    ! spherical part of the potential  &
    ! non-sph. part of the potential
    complex (kind=dp) :: rll(irmd, lmmax0d, lmmax0d), dens
    ! local variables
    real (kind=dp) :: cleb(*), drdi(irmd), visp(irmd, *), vins(irmind:irmd, lmpotd, *)
    integer :: icleb(ncleb, 4), lmsp(natypd, *), ircut(0:ipand), ipan, icell

    ! ..
    ! ---> first calculate only the spherically symmetric contribution
    real (kind=dp) :: c0ll
    complex (kind=dp), allocatable :: rsp(:), rges(:)
    integer :: lm1p, lm2p, lm3p, ir, j, i
    integer :: ircutm(0:ipand)
    ! (for all points r; if r>r_MT (or IR> IRMIN),the density has to
    ! multiplied with the shape functions...

    ! ---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)



    ! Compute spherical contribution to the torque (LM=1)
    ! Sph. potential has to be multiplied by sqrt(4 PI) !
    allocate (rges(irmd))
    allocate (rsp(irmd))

    c0ll = 1.0e0_dp/sqrt(16.0e0_dp*atan(1.0e0_dp))
    rsp = 0e0_dp
    rges = 0e0_dp

    ! cut contributions from outside the MT if recquired

    do lm1p = 1, lmmax0d
      do ir = 1, irmd
        rsp(ir) = rsp(ir) + rll(ir, lm1p, lm1p)*c0ll*(-1)*(visp(ir,nspin*(iatom-1)+2)-visp(ir,nspin*(iatom-1)+1))*0.5_dp*sqrt(16.0e0_dp*atan(1.0e0_dp))
        ! always >= 2 here
      end do
    end do

    do ir = 1, irmd
      ! --->   calculate the non spherically symmetric contribution
      if (torque_operator_onlyMT .and. (ir>ircut(1))) then
        rges(ir) = 0
      else
        rges(ir) = rsp(ir)
      end if
    end do
    if (.not. torque_operator_onlySph) then
      do j = 1, iend
        lm1p = icleb(j, 1)
        lm2p = icleb(j, 2)
        lm3p = icleb(j, 3)

        if (ipan>1 .and. lmsp(icell,lm3p)>0) then
          if (lm1p==lm2p) then
            do ir = irmin, ircut(ipan)
              rges(ir) = rges(ir) + rll(ir, lm2p, lm1p)*cleb(j)*(-1)*(vins(ir,lm3p,nspin*(iatom-1)+2)-vins(ir,lm3p,nspin*(iatom-1)+1))*0.5_dp
            end do
          else
            do ir = irmin, ircut(ipan)
              rges(ir) = rges(ir) + cleb(j)*(-1)*(vins(ir,lm3p,nspin*(iatom-1)+2)-vins(ir,lm3p,nspin*(iatom-1)+1))*0.5_dp*(rll(ir,lm2p,lm1p)+rll(ir,lm1p,lm2p))
            end do
          end if
        end if

      end do
    end if
    ! This subroutine computes a matrix that is the basis for constructing
    if (ipan==1) then
      ircutm(0) = 0
      ircutm(1) = irws
    else
      do i = 0, ipan
        ircutm(i) = ircut(i)
      end do
    end if
    ! the KKR representation of the torque operator. It is adapted from the
    call csimpk(rges(:), dens, ipan, ircutm, drdi)
    ! CALC_RHO_LL_SS subroutine, but the spin dependent part, i.e., the exhange
    deallocate (rges)
    deallocate (rsp)
    ! field, replaces the shape function in the integration.
  end subroutine calc_torq_ll_ss

end module mod_calc_torq_ll_ss
