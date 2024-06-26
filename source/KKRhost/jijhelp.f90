!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_jijhelp

  private
  public :: set_jijcalc_flags, calc_dtmatjij

contains

  !-------------------------------------------------------------------------------
  !> Summary: Set flags determining if wavefunctions and t-matrix elements for Jijs are computed
  !> Author:
  !> Category: KKRhost, physical-observables
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !-------------------------------------------------------------------------------
  subroutine set_jijcalc_flags(t_dtmatjij, natypd, natomimpd, natomimp, atomimp, iqat)

    use :: mod_types, only: t_inc, type_dtmatjijdij
    use :: mod_save_wavefun, only: t_wavefunctions

    implicit none
    integer, intent (in) :: natypd, natomimpd, natomimp, atomimp(natomimpd), iqat(natypd)
    type (type_dtmatjijdij), intent (inout) :: t_dtmatjij(natypd)

    integer :: i1, iq

    do iq = 1, natypd
      do i1 = 1, natomimp
        if (iqat(iq)==atomimp(i1)) then
          t_dtmatjij(iq)%calculate = .true.
        end if
      end do                       ! I1
    end do                         ! IQ

    ! change saveing option for wavefunctions
    t_wavefunctions%save_rll = .true.
    t_wavefunctions%save_sll = .true.
    t_wavefunctions%save_rllleft = .true.
    t_wavefunctions%save_sllleft = .true.

    ! Test
    ! write(*,*) '=========== TEST FOR TMAT-CALC============='
    ! do IQ=1,t_inc%NATYP
    ! write(*,'(A,I8,L3)') "atom",IQ, t_dtmatJij(IQ)%calculate
    ! end do

  end subroutine set_jijcalc_flags



  !-------------------------------------------------------------------------------
  !> Summary: Compute \[\Delta t]\-matrix for Jijs
  !> Author:
  !> Category: KKRhost, physical-observables, single-site
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !-------------------------------------------------------------------------------
  subroutine calc_dtmatjij(lmmax0d, lmmaxd, lmpotd, ntotd, nrmaxd, nsra, irmdnew, nspin, vins, rllleft, rll, rpan_intervall, ipan_intervall, npan_tot, ncheb, cleb, icleb, iend, &
    ncleb, rnew, dtmat)
    ! subroutine
    ! calc_dtmatJij(NTOTD,NRMAXD,NSRA,IRMDNEW,NSPIN,VINS,RLLLEFT,RLL,RPAN_INTERVALL,IPAN_INTERVALL,NPAN_TOT,NCHEB,CLEB,ICLEB,IEND,NCLEB,RNEW,dtmat)

    use :: mod_datatypes, only: dp
    use :: mod_vllmat, only: vllmat
    use :: mod_intcheb_cell, only: intcheb_cell
    use :: mod_constants, only: czero, cone
    implicit none

    integer, intent (in) :: lmmax0d, lmmaxd, lmpotd, nsra, irmdnew, nrmaxd, nspin, iend, ncleb, ntotd ! integer arguments that only define array sizes
    ! integer, intent(in) :: NSRA,IRMDNEW,NRMAXD,NSPIN,IEND,NCLEB,NTOTD  !integer arguments that only define array sizes
    integer, intent (in) :: npan_tot, ncheb, ipan_intervall(0:ntotd), icleb(ncleb, 4) ! integer arguments
    real (kind=dp), intent (in) :: rpan_intervall(0:ntotd), vins(irmdnew, lmpotd, nspin), cleb(ncleb), rnew(nrmaxd)
    complex (kind=dp), intent (in) :: rll(nsra*lmmaxd, lmmaxd, irmdnew), rllleft(nsra*lmmaxd, lmmaxd, irmdnew)
    complex (kind=dp), intent (out) :: dtmat(lmmaxd, lmmaxd, 3)

    ! ..locals..
    integer :: ispin1, ispin2, jspin, ir, ishift1, ishift2, lm1, lm2
    real (kind=dp) :: bins(irmdnew, lmpotd, 1)
    complex (kind=dp) :: bnspll0(lmmax0d, lmmax0d, irmdnew), sigma(2, 2, 3), vnspll0(lmmaxd, lmmaxd), pnsil(lmmaxd, lmmaxd), pnsir(lmmaxd, lmmaxd), cmattmp(lmmaxd, lmmaxd), &
      wr(lmmaxd, lmmaxd, irmdnew), ctmp(irmdnew)


    if (nspin/=2) stop 'calc_dtmatJij: NSPIN must be 2'

    ! BINS(:,:,1) = (VINS(:,:,2)-VINS(:,:,1))/2 !Bxc-field  !CLEANUP: how to choose the sign here????
    bins(:, :, 1) = (vins(:,:,1)-vins(:,:,2))/2 ! Bxc-field


    ! convert B_L into B_LL' by using the Gaunt coefficients
    bnspll0 = czero
    call vllmat(1, nrmaxd, irmdnew, lmmax0d, lmmax0d, bnspll0, bins, lmpotd, cleb, icleb, iend, 1, 0.0_dp, rnew, 0, ncleb)

    ! get the pauli spin matrices
    call calc_sigma(sigma)         ! use this to perform automatically infinitesimal rotations


    ! loop over sigma_{x,y,z}
    do jspin = 1, 2                ! x,y,z
      do ir = 1, irmdnew

        ! construct sigma*B_LL'(r) = VNSPLL0
        vnspll0 = czero
        do ispin1 = 1, nspin
          ishift1 = (ispin1-1)*lmmax0d
          do ispin2 = 1, nspin
            ishift2 = (ispin2-1)*lmmax0d
            vnspll0(ishift2+1:ishift2+lmmax0d, ishift1+1:ishift1+lmmax0d) = sigma(ispin2, ispin1, jspin)*bnspll0(:, :, ir)
          end do                   ! ispin2
        end do                     ! ispin1

        pnsir(:, :) = rll(1:lmmaxd, :, ir)
        pnsil(:, :) = rllleft(1:lmmaxd, :, ir)

        ! calculate [Rleft * VNSPLL0 *Rright](r)
        call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, cone, vnspll0, lmmaxd, pnsir, lmmaxd, czero, cmattmp, lmmaxd)

        call zgemm('T', 'N', lmmaxd, lmmaxd, lmmaxd, cone, pnsil, lmmaxd, cmattmp, lmmaxd, czero, wr(:,:,ir), lmmaxd)

        ! add small component
        if (nsra==2) then
          pnsir(:, :) = rll(lmmaxd+1:2*lmmaxd, :, ir)
          pnsil(:, :) = -rllleft(lmmaxd+1:2*lmmaxd, :, ir)

          ! calculate [Rleft_small * VNSPLL0 *Rright_small](r) and add to large component
          call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, cone, vnspll0, lmmaxd, pnsir, lmmaxd, czero, cmattmp, lmmaxd)

          call zgemm('T', 'N', lmmaxd, lmmaxd, lmmaxd, cone, pnsil, lmmaxd, cmattmp, lmmaxd, cone, wr(:,:,ir), lmmaxd)
        end if                     ! NSRA

      end do                       ! ir

      ! perform radial integration for each matrix element dtmat(LM2,LM1) = \int dr {Rleft * VNSPLL0 *Rright}(r)
      do lm1 = 1, lmmaxd
        do lm2 = 1, lmmaxd
          ctmp = wr(lm2, lm1, :)
          call intcheb_cell(ctmp, dtmat(lm2,lm1,jspin), rpan_intervall, ipan_intervall, npan_tot, ncheb, irmdnew)
        end do                     ! LM2
      end do                       ! LM1

    end do                         ! jspin

  end subroutine calc_dtmatjij


  ! subroutine calclambda(lambda, theta, phi)
  ! use mod_rotatespinframe, only: rotatematrix
  ! implicit none
  ! complex (kind=dp) :: lambda(2, 2, 3)
  ! real (kind=dp) :: theta, phi
  ! complex (kind=dp) :: sigmatemp(2, 2, 3), sigma(2, 2, 3)
  ! integer :: ispin

  ! call calc_sigma(sigma)
  ! sigmatemp = sigma
  ! do ispin = 1, 3
  ! call rotatematrix(sigmatemp(:,:,ispin), theta, phi, 1, 1)
  ! lambda(:, :, ispin) = sigmatemp(:, :, ispin) ! test
  ! end do                         ! ispin
  ! end subroutine calclambda


  !-------------------------------------------------------------------------------
  !> Summary: Calculate Pauli matrices
  !> Author:
  !> Category: KKRhost, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !-------------------------------------------------------------------------------
  subroutine calc_sigma(sigma)
    use :: mod_datatypes, only: dp
    implicit none
    complex (kind=dp) :: sigma(2, 2, 3)
    integer :: verbose
    character (len=*), parameter :: conventionmode = 'kkr'

    verbose = 0

    if (conventionmode=='normal') then
      sigma(1, 1, 1) = (0.0d0, 0.0d0)
      sigma(1, 2, 1) = (1.0d0, 0.0d0)
      sigma(2, 1, 1) = (1.0d0, 0.0d0)
      sigma(2, 2, 1) = (0.0d0, 0.0d0)

      sigma(1, 1, 2) = (0.0d0, 0.0d0)
      sigma(1, 2, 2) = (0.0d0, -1.0d0)
      sigma(2, 1, 2) = (0.0d0, 1.0d0)
      sigma(2, 2, 2) = (0.0d0, 0.0d0)

      sigma(1, 1, 3) = (1.0d0, 0.0d0)
      sigma(1, 2, 3) = (0.0d0, 0.0d0)
      sigma(2, 1, 3) = (0.0d0, 0.0d0)
      sigma(2, 2, 3) = (-1.0d0, 0.0d0)
    else if (conventionmode=='kkr') then
      sigma(1, 1, 1) = (0.0d0, 0.0d0)
      sigma(1, 2, 1) = (1.0d0, 0.0d0)
      sigma(2, 1, 1) = (1.0d0, 0.0d0)
      sigma(2, 2, 1) = (0.0d0, 0.0d0)

      sigma(1, 1, 2) = (0.0d0, 0.0d0)
      sigma(1, 2, 2) = (0.0d0, 1.0d0)
      sigma(2, 1, 2) = (0.0d0, -1.0d0)
      sigma(2, 2, 2) = (0.0d0, 0.0d0)

      sigma(1, 1, 3) = (-1.0d0, 0.0d0)
      sigma(1, 2, 3) = (0.0d0, 0.0d0)
      sigma(2, 1, 3) = (0.0d0, 0.0d0)
      sigma(2, 2, 3) = (1.0d0, 0.0d0)
    else
      stop '[calc_sigma] wrong mode'
    end if

    if (verbose==1) then
      write (*, *) '#################################'
      write (*, *) 'calculation OF Pauli matricies'
      write (*, *) '#################################'
      write (*, *) 'sigma_x'
      write (*, '(4f6.2)') sigma(1, 1, 1), sigma(1, 2, 1)
      write (*, '(4f6.2)') sigma(2, 1, 1), sigma(2, 2, 1)

      write (*, *) 'sigma_y'
      write (*, '(4f6.2)') sigma(1, 1, 2), sigma(1, 2, 2)
      write (*, '(4f6.2)') sigma(2, 1, 2), sigma(2, 2, 2)

      write (*, *) 'sigma_z'
      write (*, '(4f6.2)') sigma(1, 1, 3), sigma(1, 2, 3)
      write (*, '(4f6.2)') sigma(2, 1, 3), sigma(2, 2, 3)
      write (*, *) '#################################'
    end if

  end subroutine calc_sigma

end module mod_jijhelp
