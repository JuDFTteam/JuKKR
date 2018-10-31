!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Calculation of the density for the new solver
!> Author:
!> Calculation of the density for the new solver
!------------------------------------------------------------------------------------
module mod_rhooutnew

contains

  !> Summary: Calculation of the density for the new solver
  !> Author: 
  !> Category: physical-observables, KKRhost
  !> Deprecated: False 
  !> Calculation of the density for the new solver
  subroutine rhooutnew(nsra, lmax, gmatll, ek, lmpot, df, npan_tot, ncheb, cleb, icleb, iend, irmdnew, thetasnew, ifunm, imt1, lmsp, rll, rllleft, sllleft, cden, cdenlm, cdenns, &
    rho2nsc, corbital, gflle_part, rpan_intervall, ipan_intervall, nspin)

    use :: mod_constants, only: cone,czero,pi
    use :: mod_profiling
    use :: global_variables
    use :: mod_datatypes, only: dp
    use :: mod_orbitalmoment, only: calc_orbitalmoment
    use :: mod_intcheb_cell

    implicit none

    integer, intent (in) :: nsra
    integer, intent (in) :: nspin
    integer, intent (in) :: lmax   !! Maximum l component in wave function expansion
    integer, intent (in) :: iend   !! Number of nonzero gaunt coefficients
    integer, intent (in) :: imt1
    integer, intent (in) :: ncheb  !! Number of Chebychev pannels for the new solver
    integer, intent (in) :: lmpot  !! (LPOT+1)**2
    integer, intent (in) :: irmdnew
    integer, intent (in) :: corbital
    integer, intent (in) :: npan_tot
    complex (kind=dp), intent (in) :: ek
    complex (kind=dp), intent (in) :: df
    integer, dimension (*), intent (in) :: lmsp !! 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
    integer, dimension (*), intent (in) :: ifunm
    integer, dimension (0:ntotd), intent (in) :: ipan_intervall
    integer, dimension (ncleb, 4), intent (in) :: icleb !! Pointer array
    real (kind=dp), dimension (*), intent (in) :: cleb !! GAUNT coefficients (GAUNT)
    real (kind=dp), dimension (0:ntotd), intent (in) :: rpan_intervall
    real (kind=dp), dimension (ntotd*(ncheb+1), nfund), intent (in) :: thetasnew
    complex (kind=dp), dimension (lmmaxso, lmmaxso), intent (in) :: gmatll !! GMATLL=diagonal elements of the G matrix (system) Note that SLL is not needed for calculation of density, only needed for calculation of Green function
    complex (kind=dp), dimension (nsra*lmmaxso, lmmaxso, irmdnew), intent (in) :: rll
    complex (kind=dp), dimension (nsra*lmmaxso, lmmaxso, irmdnew), intent (in) :: rllleft
    complex (kind=dp), dimension (nsra*lmmaxso, lmmaxso, irmdnew), intent (in) :: sllleft

    ! .. Output variables
    complex (kind=dp), dimension (irmdnew, nspin*(1+korbit)), intent (out) :: cdenns
    complex (kind=dp), dimension (lmmaxso, lmmaxso), intent (out) :: gflle_part
    ! lmlm-dos
    complex (kind=dp), dimension (irmdnew, 0:lmax, nspin*(1+korbit)), intent (out) :: cden
    complex (kind=dp), dimension (irmdnew, lmmaxd/(1+korbit), nspin*(1+korbit)), intent (out) :: cdenlm

    ! .. In/Out variables
    complex (kind=dp), dimension (irmdnew, lmpot, nspin*(1+korbit)), intent (out) :: rho2nsc

    ! .. Local variables
    integer :: lmsize
    integer :: ir, jspin, lm1, lm2, lm3, m1, l1, j, ifun
    integer :: i_stat, i_all
    real (kind=dp) :: c0ll
    complex (kind=dp) :: cltdf
    integer, dimension (4) :: lmshift1
    integer, dimension (4) :: lmshift2
    complex (kind=dp), dimension (lmmaxso, lmmaxso, 3) :: loperator

    ! .. Local allocatable arrays
    complex (kind=dp), dimension (:), allocatable :: cwr ! lmlm-dos
    complex (kind=dp), dimension (:, :), allocatable :: qnsi
    complex (kind=dp), dimension (:, :), allocatable :: pnsi
    complex (kind=dp), dimension (:, :, :), allocatable :: wr
    complex (kind=dp), dimension (:, :, :), allocatable :: wr1 ! LDAU
    ! .. External routines
    logical, external :: opt, test

    lmsize = lmmaxd/(1+korbit)

    allocate (wr(lmmaxso,lmmaxso,irmdnew), stat=i_stat)
    call memocc(i_stat, product(shape(wr))*kind(wr), 'WR', 'RHOOUTNEW')
    wr = czero
    allocate (cwr(irmdnew), stat=i_stat)
    call memocc(i_stat, product(shape(cwr))*kind(cwr), 'CWR', 'RHOOUTNEW')
    cwr = czero
    allocate (wr1(lmmaxso,lmmaxso,irmdnew), stat=i_stat)
    call memocc(i_stat, product(shape(wr1))*kind(wr1), 'WR1', 'RHOOUTNEW')
    wr1 = czero
    allocate (qnsi(lmmaxso,lmmaxso), stat=i_stat)
    call memocc(i_stat, product(shape(qnsi))*kind(qnsi), 'QNSI', 'RHOOUTNEW')
    qnsi = czero
    allocate (pnsi(lmmaxso,lmmaxso), stat=i_stat)
    call memocc(i_stat, product(shape(pnsi))*kind(pnsi), 'PNSI', 'RHOOUTNEW')
    pnsi = czero

    ! set LMSHIFT value which is need to construct CDEN
    if (test('NOSOC   ')) then
      lmshift1(:) = 0 
      lmshift2(:) = 0 
    else
      lmshift1(1) = 0
      lmshift1(2) = lmsize
      lmshift1(3) = 0
      lmshift1(4) = lmsize
      lmshift2(1) = 0
      lmshift2(2) = lmsize
      lmshift2(3) = lmsize
      lmshift2(4) = 0
    end if

    ! for orbital moment
    if (corbital/=0) then
      call calc_orbitalmoment(lmax, lmmaxso, loperator)
    end if

    c0ll = 1e0_dp/sqrt(4e0_dp*pi)
    cden = czero
    cdenlm = czero

    ! big component of Dirac spinor
    do ir = 1, irmdnew
      do lm1 = 1, lmmaxso
        do lm2 = 1, lmmaxso
          qnsi(lm1, lm2) = sllleft(lm1, lm2, ir)
          ! PNSI(LM1,LM2)=RLL(LM1,LM2,IR)
          pnsi(lm1, lm2) = rllleft(lm1, lm2, ir)
        end do
      end do
      ! CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,PNSI,
      ! +             LMMAXSO,GMATLL,LMMAXSO,EK,QNSI,LMMAXSO)
      call zgemm('N', 'T', lmmaxso, lmmaxso, lmmaxso, cone, pnsi, lmmaxso, gmatll, lmmaxso, ek, qnsi, lmmaxso)
      do lm1 = 1, lmmaxso
        do lm2 = 1, lmmaxso
          pnsi(lm1, lm2) = rll(lm1, lm2, ir)
        end do
      end do
      call zgemm('N', 'T', lmmaxso, lmmaxso, lmmaxso, cone, pnsi, lmmaxso, qnsi, lmmaxso, czero, wr(1,1,ir), lmmaxso)

      ! small component of Dirac spinor
      if (nsra==2) then
        do lm1 = 1, lmmaxso
          do lm2 = 1, lmmaxso
            ! QNSI(LM1,LM2)=SLLLEFT(LM1+LMMAXSO,LM2,IR)
            qnsi(lm1, lm2) = -sllleft(lm1+lmmaxso, lm2, ir)
            ! PNSI(LM1,LM2)=RLLLEFT(LM1+LMMAXSO,LM2,IR)
            pnsi(lm1, lm2) = -rllleft(lm1+lmmaxso, lm2, ir)
          end do
        end do
        ! CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,PNSI,
        ! +             LMMAXSO,GMATLL,LMMAXSO,EK,QNSI,LMMAXSO)
        call zgemm('N', 'T', lmmaxso, lmmaxso, lmmaxso, cone, pnsi, lmmaxso, gmatll, lmmaxso, ek, qnsi, lmmaxso)
        do lm1 = 1, lmmaxso
          do lm2 = 1, lmmaxso
            pnsi(lm1, lm2) = rll(lm1+lmmaxso, lm2, ir)
          end do
        end do
        call zgemm('N', 'T', lmmaxso, lmmaxso, lmmaxso, cone, pnsi, lmmaxso, qnsi, lmmaxso, cone, wr(1,1,ir), lmmaxso)
      end if ! small component

      ! For orbital moment
      if (corbital/=0) then
        call zgemm('N', 'N', lmmaxso, lmmaxso, lmmaxso, cone, loperator(1,1,corbital), lmmaxso, wr(1,1,ir), lmmaxso, czero, pnsi, lmmaxso)
        do lm1 = 1, lmmaxso
          do lm2 = 1, lmmaxso
            wr(lm1, lm2, ir) = pnsi(lm1, lm2)
          end do
        end do
      end if
      do lm1 = 1, lmmaxso
        do lm2 = 1, lmmaxso
          wr1(lm1, lm2, ir) = wr(lm1, lm2, ir)
        end do
      end do
      do lm1 = 1, lmmaxso
        do lm2 = 1, lm1 - 1
          wr1(lm1, lm2, ir) = wr1(lm1, lm2, ir) + wr1(lm2, lm1, ir)
        end do
      end do

      do jspin = 1, nspin*(1+korbit)
        do lm1 = 1, lmsize
          do lm2 = 1, lm1 - 1
            wr(lm1+lmshift1(jspin), lm2+lmshift2(jspin), ir) =                      &
              wr(lm1+lmshift1(jspin), lm2+lmshift2(jspin), ir) +                    &
              wr(lm2+lmshift1(jspin), lm1+lmshift2(jspin), ir)
          end do
        end do
      end do ! JSPIN

    end do ! IR


    ! IF lmdos or LDAU
    if (opt('lmlm-dos') .or. opt('LDA+U   ')) then ! lmlm-dos
      ! Integrate only up to muffin-tin radius.
      ! ! lmlm-dos
      gflle_part = czero           ! lmlm-dos
      do lm2 = 1, lmmaxso          ! lmlm-dos
        do lm1 = 1, lmmaxso        ! lmlm-dos
          ! For integration up to MT radius do this:                           !
          ! lmlm-dos
          ! CWR(1:IMT1) = WR(LM1,LM2,1:IMT1)                                   !
          ! lmlm-dos
          ! CWR(IMT1+1:IRMDNEW) = CZERO                                        !
          ! lmlm-dos
          ! CALL INTCHEB_CELL(CWR,GFLLE_PART(LM1,LM2),RPAN_INTERVALL,&         !
          ! lmlm-dos
          ! IPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW)                         !
          ! lmlm-dos
          ! For full cell integration replace loop content with this:          !
          ! lmlm-dos
          cwr(1:irmdnew) = wr1(lm1, lm2, 1:irmdnew) ! lmlm-dos
          ! If LDAU, integrate only up to MT
          do ir = imt1 + 1, irmdnew
            if (opt('LDA+U   ')) then
              cwr(ir) = czero      ! LDAU
            else
              cwr(ir) = cwr(ir)*thetasnew(ir, 1)*c0ll ! lmlm-dos
            end if
          end do
          call intcheb_cell(cwr, gflle_part(lm1,lm2), rpan_intervall, ipan_intervall, npan_tot, ncheb, irmdnew)
        end do
      end do
    end if                         ! OPT('lmlm-dos').OR.OPT('LDA+U   ')


    ! DO IR = 1,IRMDNEW
    ! DO JSPIN = 1,4
    ! DO LM1 = 1,lmsize
    ! DO LM2 = 1,LM1-1
    ! WR(LM1+LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),IR)=
    ! +           WR(LM1+LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),IR)+
    ! +           WR(LM2+LMSHIFT1(JSPIN),LM1+LMSHIFT2(JSPIN),IR)
    ! ENDDO
    ! ENDDO
    ! ENDDO ! JSPIN
    ! ENDDO !IR

    ! First calculate the spherical symmetric contribution

    do l1 = 0, lmax
      do m1 = -l1, l1
        lm1 = l1*(l1+1) + m1 + 1
        do ir = 1, irmdnew
          do jspin = 1, nspin*(1+korbit)
            cden(ir, l1, jspin) = cden(ir, l1, jspin) + wr(lm1+lmshift1(jspin), lm1+lmshift2(jspin), ir)
            cdenlm(ir, lm1, jspin) = wr(lm1+lmshift1(jspin), lm1+lmshift2(jspin), ir)
          end do                   ! JPSIN
        end do                     ! IR
      end do                       ! M1

      do jspin = 1, nspin*(1+korbit)
        do ir = 1, irmdnew
          rho2nsc(ir, 1, jspin) = rho2nsc(ir, 1, jspin) + c0ll*(cden(ir,l1,jspin)*df)
        end do                     ! IR

        do ir = imt1 + 1, irmdnew
          cden(ir, l1, jspin) = cden(ir, l1, jspin)*thetasnew(ir, 1)*c0ll
          do m1 = -l1, l1
            lm1 = l1*(l1+1) + m1 + 1
            cdenlm(ir, lm1, jspin) = cdenlm(ir, lm1, jspin)*thetasnew(ir, 1)*c0ll
          end do                   ! M1
        end do                     ! IR
      end do                       ! JSPIN
    end do                         ! L1

    ! Then the non-spherical part

    cdenns = czero
    do j = 1, iend
      lm1 = icleb(j, 1)
      lm2 = icleb(j, 2)
      lm3 = icleb(j, 3)
      cltdf = df*cleb(j)
      do jspin = 1, nspin*(1+korbit)
        do ir = 1, irmdnew
          rho2nsc(ir, lm3, jspin) = rho2nsc(ir, lm3, jspin) + (cltdf*wr(lm1+lmshift1(jspin),lm2+lmshift2(jspin),ir))
        end do

        if (lmsp(lm3)>0) then
          ifun = ifunm(lm3)
          do ir = imt1 + 1, irmdnew
            cdenns(ir, jspin) = cdenns(ir, jspin) + cleb(j)*wr(lm1+lmshift1(jspin), lm2+lmshift2(jspin), ir)*thetasnew(ir, ifun)
          end do
        end if
      end do                       ! JSPIN
    end do                         ! J


    i_all = -product(shape(wr))*kind(wr)
    deallocate (wr, stat=i_stat)
    call memocc(i_stat, i_all, 'WR', 'RHOOUTNEW')
    i_all = -product(shape(wr1))*kind(wr1)
    deallocate (wr1, stat=i_stat)
    call memocc(i_stat, i_all, 'WR1', 'RHOOUTNEW')
    i_all = -product(shape(cwr))*kind(cwr)
    deallocate (cwr, stat=i_stat)
    call memocc(i_stat, i_all, 'CWR', 'RHOOUTNEW')
    i_all = -product(shape(qnsi))*kind(qnsi)
    deallocate (qnsi, stat=i_stat)
    call memocc(i_stat, i_all, 'QNSI', 'RHOOUTNEW')
    i_all = -product(shape(pnsi))*kind(pnsi)
    deallocate (pnsi, stat=i_stat)
    call memocc(i_stat, i_all, 'PNSI', 'RHOOUTNEW')

  end subroutine rhooutnew

end module mod_rhooutnew
