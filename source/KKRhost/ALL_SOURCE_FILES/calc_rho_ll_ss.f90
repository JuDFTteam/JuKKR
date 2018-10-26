!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_calc_rho_ll_ss

  private
  public :: calc_rho_ll_ss

contains

  !-------------------------------------------------------------------------------
  !> Summary: Rho matrix for PKKprime code
  !> Author: 
  !> Category: KKRhost, physical-observables
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !> 
  !> Computes rho matrix (see eq. 2.111 ff. in PhD thesis of Bernd Zimmermann)
  !-------------------------------------------------------------------------------
  subroutine calc_rho_ll_ss(lmmax, rll, ircut, ipan, icell, thetas, cleb, icleb, iend, ifunm, lmsp, irws, drdi, dens)

    use :: mod_datatypes, only: dp
    use :: global_variables, only: irmd, irid, nfund, ncleb, natypd, lmpotd, ipand
    use :: mod_csimpk, only: csimpk
    implicit none

    ! .. Array Arguments ..
    ! non-sph. eigen states of single pot
    integer :: iend, lmmax, irws   ! derivative dr/di

    ! local variables
    ! ,DENS(:,:,:)
    complex (kind=dp) :: rll(irmd, lmmax, lmmax), dens
    real (kind=dp) :: cleb(*), thetas(irid, nfund, *), drdi(irmd) ! RGES_W(:,:,:,:),
    ! &
    integer :: icleb(ncleb, 4), ifunm(natypd, lmpotd), lmsp(natypd, *), ircut(0:ipand), ipan, icell, ifun
    ! DENS_GESAMT(:,:), &
    ! DENS_GESAMT_I1(:,:,:)
    real (kind=dp) :: c0ll
    complex (kind=dp) :: clt
    complex (kind=dp), allocatable :: rsp(:), rges(:) ! ..
    ! ---> first calculate only the spherically symmetric contribution
    ! (for all points r; if r>r_MT (or IR> IRMIN),the density has to
    ! multiplied with the shape functions...
    integer :: lm1p, lm2p, lm3p, ir, j, i
    integer :: ircutm(0:ipand)

    ! ---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)

    ! WRITE(6,*) "In rho ll"

    allocate (rges(irmd))
    allocate (rsp(irmd))

    ! WRITE(56,"((I5),(2e17.9))") IR,THETAS(IR-IRCUT(1),1,ICELL)
    c0ll = 1.0e0_dp/sqrt(16.0e0_dp*atan(1.0e0_dp))
    rsp = 0e0_dp
    rges = 0e0_dp

    do lm1p = 1, lmmax
      do ir = 1, irmd
        rsp(ir) = rsp(ir) + rll(ir, lm1p, lm1p)
      end do
    end do
    ! STOP " "
    do ir = 1, ircut(ipan)
      rges(ir) = rsp(ir)
    end do
    ! WRITE(6,*) "IRCUT(1)",IRCUT(1)
    if (ipan>1) then
      do ir = ircut(1) + 1, ircut(ipan)
        ! WRITE(6,*) "IRCUT(IPAN)",IRCUT(IPAN)
        rges(ir) = rsp(ir)*c0ll*thetas(ir-ircut(1), 1, icell)
      end do
    end if
    ! WRITE(6,*) "IRCUT(IPAN)-IRMIND",IRCUT(IPAN)-IRMIND
    ! WRITE(6,*) "IRMIND",IRMIND

    ! ---> calculate the non spherically symmetric contribution

    ! WRITE(156,*) "IFUN",IFUN
    do j = 1, iend
      lm1p = icleb(j, 1)
      lm2p = icleb(j, 2)
      lm3p = icleb(j, 3)
      clt = cleb(j)

      if (ipan>1 .and. lmsp(icell,lm3p)>0) then
        ifun = ifunm(icell, lm3p)

        if (lm1p==lm2p) then
          do ir = ircut(1) + 1, ircut(ipan)
            rges(ir) = rges(ir) + rll(ir, lm2p, lm1p)*cleb(j)*thetas(ir-ircut(1), ifun, icell)
          end do
        else
          do ir = ircut(1) + 1, ircut(ipan)
            rges(ir) = rges(ir) + cleb(j)*thetas(ir-ircut(1), ifun, icell)*(rll(ir,lm2p,lm1p)+rll(ir,lm1p,lm2p))
          end do
        end if

      end if

    end do

    if (ipan==1) then
      ircutm(0) = 0
      ircutm(1) = irws
    else
      do i = 0, ipan
        ircutm(i) = ircut(i)
      end do
    end if

    call csimpk(rges(:), dens, ipan, ircutm, drdi)
    ! set to 1 if NEWSOSOL under RUNOPT, otherwise 0
    deallocate (rges)
    deallocate (rsp)
    ! SET ACCORDING TO lmax VALUE OF INPUTCARD
  end subroutine calc_rho_ll_ss

end module mod_calc_rho_ll_ss
