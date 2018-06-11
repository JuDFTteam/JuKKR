!**********************************************************************
subroutine mssinit(ncpa, icpastart, tsst, msst, mssq, trefll, drotq, refpot, &
  iqat, itoq, noq, conc, kmrot, natyp, naez, lmmaxd) ! nrefd was taken out of calling list 1.2.2012

  use :: mod_mympi, only: myrank, master
      Use mod_datatypes, Only: dp
  implicit none
  include 'inc.p' ! Included  1.2.2012

!.. Local variables
  complex (kind=dp) :: czero, cone
  parameter (czero=(0.0d0,0.0d0))
  parameter (cone=(1.0d0,0.0d0))
!..
!.. External Subroutines ..
  integer :: kmrot, natyp, naez, lmmaxd, ncpa, icpastart
  integer :: iqat(natypd), itoq(natypd, naezd)
  integer :: refpot(naezd), noq(naezd)
  real (kind=dp) :: conc(natypd)
  complex (kind=dp) :: tsst(lmmaxd, lmmaxd, natypd), trefll(lmmaxd, lmmaxd, nrefd &
    )
  complex (kind=dp) :: msst(lmmaxd, lmmaxd, natypd)
  complex (kind=dp) :: mssq(lmmaxd, lmmaxd, naezd)
  complex (kind=dp) :: drotq(lmmaxd, lmmaxd, naezd)
! ======================================================================

  integer :: it, iq, rf, j, io, info, lm1, lm2, lp, ld, lmp, lmd
  integer :: ipvt(lmmaxd)
  complex (kind=dp) :: zc
  complex (kind=dp) :: w1(lmmaxd, lmmaxd)
  complex (kind=dp) :: w2(lmmaxd, lmmaxd)
!   --> set up the Delta_t^-1 matrix (MSST) in the LOCAL frame

  logical :: test, opt
  external :: test, opt







  do it = 1, natyp

    do j = 1, lmmaxd
      call zcopy(lmmaxd, tsst(1,j,it), 1, msst(1,j,it), 1)
    end do

    iq = iqat(it)
    rf = refpot(iq)
! ---> determine Delta_t = t(sys) - tmat(ref) = TSST - TREFLL
!      in local frame
    if (kmrot/=0) then
      call rotate(trefll(1,1,rf), 'G->L', w1, lmmaxd, drotq(1,1,iq), lmmaxd)
    else
      do j = 1, lmmaxd
        call zcopy(lmmaxd, trefll(1,j,rf), 1, w1(1,j), 1)
      end do
    end if




!  --> inversion
    do lm1 = 1, lmmaxd
      do lm2 = 1, lmmaxd

        msst(lm2, lm1, it) = (msst(lm2,lm1,it)-w1(lm2,lm1))
!         CALL ZGETRI(LMMAXD,MSST(1,1,IT),LMMAXD,IPVT,W1,
      end do
    end do
!     &               LMMAXD*LMMAXD,INFO)
!( TEST('testgmat') ) THEN
!( OPT('VIRATOMS') ) THEN
    if (.not. opt('VIRATOMS')) then
      if (.not. test('testgmat')) then
        call zgetrf(lmmaxd, lmmaxd, msst(1,1,it), lmmaxd, ipvt, info)


        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            if (lm1==lm2) then
              w2(lm1, lm2) = cone
            else
              w2(lm1, lm2) = czero
            end if
          end do
        end do
        call zgetrs('N', lmmaxd, lmmaxd, msst(1,1,it), lmmaxd, ipvt, w2, &
          lmmaxd, info)
        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            msst(lm1, lm2, it) = w2(lm1, lm2)
          end do
        end do
      end if 
    end if ! ======================================================================

  end do
! ---> determine tmat(sys) - tmat(ref) = MSSQ - TREFLL

!      because TSST is calculated in the LOCAL frame, if KMROT<>0
!      it needs to be rotated prior to set the Delta matrix

! ---> set up the effective (on-site) Delta_t- and Delta_m-matrices
!      using the Average T-matrix Approximation

!   ICPASTART=1:
!       m(IQ) = t(ATA) = SUM(it)  c(it) * t(it)
!       m(IQ) = t(ATA) - t_ref = Delta_t(ATA)
!       m(IQ) = (Delta_t(ATA))^(-1)

!   ICPASTART=2:
!       m(IQ) = (Delta_t(ATA))^(-1) for l = 2
!       m(IQ) = SUM(it) c(it)*m(it) with m(it)=t(it)^(-1)

!       mssq(IQ)  refer to the GLOBAL frame
!       tsst(IT),msst(IT)  refer to the LOCAL  frame


! ----------------------------------------------------------------------



!         write(*,*) 'test fivos mssinit IO,IQ,IT',IO,IQ,IT ! test fivos

  call cinit(lmmaxd*lmmaxd*naez, mssq)
! ---> rotate MSSQ from the LOCAL to the GLOBAL frame if necessary
  do iq = 1, naez

    do io = 1, noq(iq)
      it = itoq(io, iq)

      zc = conc(it)
      do j = 1, lmmaxd
        call zaxpy(lmmaxd, zc, tsst(1,j,it), 1, mssq(1,j,iq), 1)
      end do
    end do



    if (kmrot/=0) then
      call rotate(mssq(1,1,iq), 'L->G', w1, lmmaxd, drotq(1,1,iq), lmmaxd)
      do j = 1, lmmaxd
        call zcopy(lmmaxd, w1(1,j), 1, mssq(1,j,iq), 1)
      end do
    end if
! ---> determine Delta_t = t(sys) - tmat(ref) = TSSQ - TREFLL
!      in the GLOBAL frame

!           write(*,*) 'RF',RF



    rf = refpot(iq)


    do lm1 = 1, lmmaxd
      do lm2 = 1, lmmaxd

        mssq(lm2, lm1, iq) = (mssq(lm2,lm1,iq)-trefll(lm2,lm1,rf))
! ----------------------------------------------------------------------
      end do
    end do
!    store the Delta_t matrix
    if (test('tmat    ')) then
      write (1337, *) 'IQ,IT,RF', iq, it, rf
      write (1337, *) 'DELTA_TMATLL (', iq, ' )'
      call cmatstr(' ', 1, mssq(1,1,iq), lmmaxd, lmmaxd, 2*krel+1, 2*krel+1, &
        0, 1d-8, 6)
      write (1337, *)
    end if
! ----------------------------------------------------------------------
! fswrt
  end do
! fswrt
! fswrt
! fswrt
  if (opt('FERMIOUT') .and. myrank==master) then ! fswrt
    write (6801, '(A)') 'TMATLL(ie):' ! fswrt
    do iq = 1, naez ! fswrt
      do lm2 = 1, lmmaxd ! fswrt
        do lm1 = 1, lmmaxd ! fswrt
          write (6801, '(2ES25.16)') mssq(lm1, lm2, iq) ! fswrt
        end do ! ----------------------------------------------------------------------
      end do 
    end do !    MSSQ is now the Delta_t matrix in the GLOBAL frame
  end if !    below, we determine (Delta_t)^(-1) in the GLOBAL frame

! ----------------------------------------------------------------------

!=======================================================================

! ---> loop over all atoms in unit cell, get Delta_t^(-1) = MSSQ


! ---> inversion


  do iq = 1, naez
!        CALL ZGETRI(LMMAXD,MSSQ(1,1,IQ),LMMAXD,IPVT,W1,
!     &              LMMAXD*LMMAXD,INFO)
!( .not. TEST('testgmat') ) THEN
!( .not. OPT('VIRATOMS') ) THEN
    if (.not. opt('VIRATOMS')) then
      if (.not. test('testgmat')) then
        call zgetrf(lmmaxd, lmmaxd, mssq(1,1,iq), lmmaxd, ipvt, info)


        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            if (lm1==lm2) then
              w2(lm1, lm2) = cone
            else
              w2(lm1, lm2) = czero
            end if
          end do
        end do
        call zgetrs('N', lmmaxd, lmmaxd, mssq(1,1,iq), lmmaxd, ipvt, w2, &
          lmmaxd, info)
        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            mssq(lm1, lm2, iq) = w2(lm1, lm2)
          end do
        end do
      end if ! IQ = 1,NAEZ
    end if !            stop


  end do !============================================================IQ = 1,NAEZ
!----------------------------------------------------------------------
! s-, p-, and f-terms:    m(ata) = sum(q) c(q) * m(q)         >>>  AKAI

! ------------------------------------------ s,p blocks
  if ((ncpa/=0) .and. (icpastart==2)) then



    lp = 1
    ld = 2
    lmp = (krel+1)*(lp+1)**2
    lmd = (krel+1)*(ld+1)**2 + 1
    do iq = 1, naez
! ------------------------------------------ f block
      do lm1 = 1, lmp
        do lm2 = 1, lmp
          mssq(lm1, lm2, iq) = czero

          do io = 1, noq(iq)
            it = itoq(io, iq)
            mssq(lm1, lm2, iq) = mssq(lm1, lm2, iq) + &
              conc(it)*msst(lm1, lm2, it)
          end do

        end do
      end do

! IQ=1,NAEZ
      do lm1 = lmd, lmmaxd
        do lm2 = lmd, lmmaxd
          mssq(lm1, lm2, iq) = czero

          do io = 1, noq(iq)
            it = itoq(io, iq)
            mssq(lm1, lm2, iq) = mssq(lm1, lm2, iq) + &
              conc(it)*msst(lm1, lm2, it)
          end do
! ICPASTART.EQ.2
        end do
      end do

    end do !**********************************************************************
! nrefd was taken out of calling list 1.2.2012
  end if 
! Included  1.2.2012
  return
end subroutine
