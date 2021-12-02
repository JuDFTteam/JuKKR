module mod_rhooutnew

contains

!-------------------------------------------------------------------------
!> Summary: Calculation of valence charge density, new solver
!> Category: physical-observables, KKRimp
!>
!> @warning Changed behavior of intcheb call for LDA+U run since this was
!> probably included in ir loop by accident! @endwarning
!-------------------------------------------------------------------------
subroutine rhooutnew(gauntcoeff, df, gmatin, ek, cellnew, wavefunction, rho2nsc, &
                     nsra, lmaxd, lmaxatom, lmmaxatom, lmsize, lmsize2, lmpotd, irmd,&
                     ispin, nspinden, imt1, cden, cdenlm, cdenns, shapefun, corbital, &
                     gflle_part)                                              ! lda+u

  use mod_datatypes, only: dp
  use mod_constants, only: czero, cone, pi
  use type_gauntcoeff, only: gauntcoeff_type
  use type_cellnew, only: cell_typenew
  use type_wavefunction, only: wavefunction_type
  use mod_physic_params, only: cvlight
  use mod_mathtools, only: transpose
  use mod_config, only: config_testflag
  use type_shapefun, only: shapefun_type
  use mod_orbitalmoment, only: calc_orbitalmoment
  use mod_intcheb_cell, only: intcheb_cell   ! lda+u
  use mod_timing, only: timing_start, timing_pause, timing_stop

  implicit none
  type(gauntcoeff_type) :: gauntcoeff
  type(cell_typenew) :: cellnew
  type(wavefunction_type) :: wavefunction
  integer :: lmaxd, lmaxatom
  integer :: lmmaxatom
  integer :: lmsize, lmsize2
  integer :: lmpotd
  integer :: irmd
! ..
! .. Scalar Arguments ..
  complex (kind=dp) :: df, ek
  integer :: imt1, nsra
  type(shapefun_type),intent(in) :: shapefun
  integer :: corbital
! ..
! .. Array Arguments ..
  complex (kind=dp) :: cden(irmd,0:lmaxd,nspinden)
  complex (kind=dp) :: cdenns(irmd,nspinden)
  complex (kind=dp) :: gmat(lmsize,lmsize)
  complex (kind=dp) :: gmatin(lmsize,lmsize)
  complex (kind=dp) :: qnsi(lmsize,lmsize), rlltemp(lmsize,lmsize)
  complex (kind=dp) :: cdenlm(irmd,lmmaxatom,nspinden) ! lm-dos
  complex (kind=dp) :: rho2nsc(irmd,lmpotd,nspinden)
  complex (kind=dp) :: gflle_part(lmsize,lmsize) ! lda+u
! ..
! .. Local Scalars ..
  complex (kind=dp) :: cltdf, alpha
  real (kind=dp) :: c0ll
  integer :: ifun, ir, j, l1, lm1, lm2, lm3, m1
  integer :: ispin, jspin, ierr
! ..
! .. Local Arrays ..
  complex (kind=dp), allocatable :: wr(:,:,:)
  complex (kind=dp), allocatable :: cwr(:) ! lda+u
  integer :: spinindex1(4) !=(/1,2,1,2 /)
  integer :: spinindex2(4) !=(/1,2,2,1 /)
  integer :: lmshift1(4)
  integer :: lmshift2(4)
  integer :: nspinstart, nspinstop, nspinden
  complex (kind=dp) :: loperator(lmsize,lmsize,3)
! ..
! .. External Subroutines ..
  external :: zgemm
! ..
! .. Intrinsic Functions ..
  intrinsic :: sqrt

  allocate ( wr(lmsize,lmsize,irmd), cwr(irmd), stat=ierr ) ! cwr for lda+u
  if (ierr/=0) stop 'Error allocating wr, cwr in rhooutnew'

  gmat = gmatin

  if (nspinden==4) then
    nspinstart = 1
    nspinstop = nspinden
    spinindex1 = (/1,2,1,2 /)
    spinindex2 = (/1,2,2,1 /)
    lmshift1 = lmmaxatom*(spinindex1-1) 
    lmshift2 = lmmaxatom*(spinindex2-1) 
  else
    nspinstart = ispin
    nspinstop = ispin
    spinindex1 =(/1,1,0,0 /)
    spinindex2 =(/1,1,0,0 /)
    lmshift1 = lmmaxatom*(spinindex1-1)
    lmshift2 = lmmaxatom*(spinindex2-1) 
  end if
     
  if (corbital/=0) then
     call calc_orbitalmoment(lmaxatom,loperator)
  end if

  c0ll = 1.0_dp/sqrt(4.0_dp*pi)

  !
  !---> initialize array for complex charge density
  !
  do jspin = nspinstart,nspinstop
     cden(:,:,jspin) = czero
     cdenlm(:,:,jspin) = czero
  end do !jspin

  !------------------------------------------------------------------
  !
  !---> set up array ek*qns(lm1,lm2) + { gmat(lm3,lm2)*pns(lm1,lm3) }
  !                                      summed over lm3
  !---> set up of wr(lm1,lm2) = { pns(lm1,lm3)*qns(lm2,lm3) }
  !                                               summed over lm3
  do ir = 1, irmd
    ! ########################################################3
    !
    ! WATCH OUT CHECK IF A FACTOR OF M_0 needs to be added into the Greensfunction
    !
    ! #########################################################3
    if (allocated(wavefunction%sllleft)) then
      qnsi(:,:) = wavefunction%sllleft(1:lmsize,1:lmsize,ir,1)
    else
      qnsi(:,:) = wavefunction%sll(1:lmsize,1:lmsize,ir,1)
    end if

    if (allocated(wavefunction%sllleft)) then
      rlltemp = wavefunction%rllleft(1:lmsize,1:lmsize,ir,1)
    else
      rlltemp = wavefunction%rll(1:lmsize,1:lmsize,ir,1)
    end if

    ! this is the prefactor for the gmatll*rllleft term in the first zgemm
    ! if the onsite densit is calculated alone we set this to zero
    alpha = cone
    if (config_testflag('calc_onsite_only')) alpha = czero

    ! changed the second mode to transpose - bauer
    call zgemm('n','t',lmsize,lmsize,lmsize, alpha,rlltemp, lmsize,gmat,lmsize,ek,qnsi,lmsize)

    rlltemp(:,:) = wavefunction%rll(1:lmsize,1:lmsize,ir,1)

    call zgemm('n','t',lmsize,lmsize,lmsize,cone,rlltemp, lmsize,qnsi,lmsize,czero,wr(1,1,ir),lmsize)


    if (nsra.eq.2 .and. (.not. config_testflag('nosmallcomp')) ) then

      if (allocated(wavefunction%sllleft)) then
         qnsi(:,:) = -wavefunction%sllleft(lmsize+1:2*lmsize,1:lmsize,ir,1) ! attention to the
                                                                            ! additional minus sign
         ! ##########################################################################################
         ! Drittler assumes that for the left solution, is given by the right solution with an
         ! additional minus sign. This minus sign is contained inside the equations to calculate
         ! the electronic density. While calculating the left solution, the minus sign is already 
         ! included in the left solution. To make calculations consistant a factor of -1 is included
         ! which cancels out by the routines of Drittler
         ! ##########################################################################################
      else
         qnsi(:,:) = wavefunction%sll(lmsize+1:2*lmsize,1:lmsize,ir,1)
      end if
      if (allocated(wavefunction%rllleft)) then
         rlltemp = -wavefunction%rllleft(lmsize+1:2*lmsize,1:lmsize,ir,1) ! attention to the
                                                                          ! additional minus sign
         ! ##########################################################################################
         ! Drittler assumes that for the left solution, is given by the right solution with an
         ! additional minus sign. This minus sign is contained inside the equations to calculate
         ! the electronic density. While calculating the left solution, the minus sign is already 
         ! included in the left solution. To make calculations consistant a factor of -1 is included
         ! which cancels out by the routines of Drittler
         ! ##########################################################################################
      else
         rlltemp = wavefunction%rll(lmsize+1:2*lmsize,1:lmsize,ir,1)
      end if

      ! changed the second mode to transpose - bauer
      call zgemm('n','t',lmsize,lmsize,lmsize,gmat,rlltemp, &
                 lmsize,gmat,lmsize,ek,qnsi,lmsize)

      rlltemp = wavefunction%rll(lmsize+1:2*lmsize,1:lmsize,ir,1)!/cvlight

      call zgemm('n','t',lmsize,lmsize,lmsize,cone,rlltemp, &
                  lmsize,qnsi,lmsize,cone,wr(1,1,ir),lmsize)

    end if

    if (corbital/=0) then
      call zgemm('n','n',lmsize,lmsize,lmsize,cone,loperator(:,:,corbital), &
                 lmsize,wr(:,:,ir),lmsize,czero,rlltemp,lmsize)
      wr(:,:,ir) = rlltemp
    end if

  end do !ir


  ! Phivos lda+u: Place here r-integration of wr(lms1,lms2,ir) to obtain cnll(lms1,lms2).    ! lda+u
  ! Integrate only up to muffin-tin radius.                                                  ! lda+u
  ! @note This was included in ir loop above which is probably not what was inteded. @endnote
  gflle_part(:,:) = czero                                                                    ! lda+u
  do lm2 = 1,lmsize                                                                          ! lda+u
    do lm1 = 1,lmsize                                                                        ! lda+u
      cwr(1:imt1) = wr(lm1,lm2,1:imt1)                                                       ! lda+u
      cwr(imt1+1:irmd) = czero                                                               ! lda+u
      call intcheb_cell(cwr,gflle_part(lm1,lm2), cellnew%rpan_intervall, cellnew%ipan_intervall, cellnew%npan_tot, cellnew%ncheb, cellnew%nrmaxnew)                           ! lda+u
    enddo                                                                                    ! lda+u
  enddo                                                                                      ! lda+u

  ! Change by Phivos 12.6.2012: this part was within previous IR-loop, now moved here        ! lda+u
  ! in order to calculate the complex density just before this averaging.                    ! lda+u
  ! @note this is a new ir loop because the above intcheb call needs to be done
  ! before this summation @endnote
  do ir = 1, irmd
    do jspin = nspinstart,nspinstop
      do lm1 = 1,lmmaxatom
        do lm2 = 1,lm1 - 1
          wr(lm1+lmshift1(jspin),lm2+lmshift2(jspin),ir) = wr(lm1+lmshift1(jspin),lm2+lmshift2(jspin),ir) + wr(lm2+lmshift1(jspin),lm1+lmshift2(jspin),ir)
        enddo
      enddo
    end do !jspin
  end do !ir

  !
  !---> first calculate only the spherically symmetric contribution
  !
  do l1 = 0,lmaxatom
    do m1 = -l1,l1
      lm1 = l1* (l1+1) + m1 + 1
      do ir = 1,irmd
        !
        !---> fill array for complex density of states
        !
        do jspin = nspinstart,nspinstop
          cden(ir,l1,jspin) = cden(ir,l1,jspin) + wr(lm1+lmshift1(jspin),lm1+lmshift2(jspin),ir)
          cdenlm(ir,lm1,jspin) = wr(lm1+lmshift1(jspin),lm1+lmshift2(jspin),ir) ! lm-dos
        end do
      end do !ir
    end do !m1
    !
    !---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
    !
    do jspin = nspinstart,nspinstop
      do ir = 1,irmd
        rho2nsc(ir,1,jspin) = rho2nsc(ir,1,jspin) + c0ll*(cden(ir,l1,jspin)*df)
      end do
      do ir = imt1 + 1,irmd
        cden(ir,l1,jspin) = cden(ir,l1,jspin)*cellnew%shapefun(ir,1)*c0ll
        do m1 = -l1,l1                                                            ! lm-dos
          lm1 = l1* (l1+1) + m1 + 1                                               ! lm-dos
          cdenlm(ir,lm1,jspin) = cdenlm(ir,lm1,jspin)*cellnew%shapefun(ir,1)*c0ll ! lm-dos
        enddo                                                                     ! lm-dos
      end do
    end do
  end do ! l1 = 0,lmaxatom

  do jspin = nspinstart,nspinstop
    cdenns(:,jspin) = 0.0_dp
  end do

  do j = 1,gauntcoeff%iend
    lm1 = gauntcoeff%icleb(j,1)
    lm2 = gauntcoeff%icleb(j,2)
    lm3 = gauntcoeff%icleb(j,3)
    cltdf = df*gauntcoeff%cleb(j,1)
    !
    !---> calculate the non spherically symmetric contribution
    !
    do jspin = nspinstart,nspinstop
      do ir = 1,irmd
        rho2nsc(ir,lm3,jspin) = rho2nsc(ir,lm3,jspin) + (cltdf*wr(lm1+lmshift1(jspin),lm2+lmshift2(jspin),ir))
      end do
      if (shapefun%lmused(lm3)==1) then
        ifun = shapefun%lm2index(lm3) !ifunm(lm3)
        do ir = imt1 + 1,irmd
          cdenns(ir,jspin) = cdenns(ir,jspin) + gauntcoeff%cleb(j,1)*wr(lm1+lmshift1(jspin),lm2+lmshift2(jspin),ir)*cellnew%shapefun(ir,ifun)
        end do
      end if
    end do

  end do !j

  deallocate(wr, cwr, stat=ierr)
  if (ierr/=0) stop 'Error deallocating wr, cwr in rhooutnew'

  end subroutine rhooutnew

end module mod_rhooutnew
