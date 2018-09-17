module mod_rhoout

contains

  ! Added IRMIN,IRMAX 1.7.2014
  subroutine rhoout(cden, df, gmat, ek, pns, qns, rho2ns, thetas, ifunm, ipan1, imt1, irmin, irmax, lmsp, cdenns, nsra, cleb, icleb, iend, cdenlm, cwr) ! lm-dos
    ! -----------------------------------------------------------------------

    ! calculates the charge density from r(irmin) to r(irc)
    ! in case of a non spherical input potential .

    ! fills the array cden for the complex density of states

    ! attention : the gaunt coeffients are stored in index array
    ! (see subroutine gaunt)

    ! the structured part of the greens-function (gmat) is symmetric in
    ! its lm-indices , therefore only one half of the matrix is
    ! calculated in the subroutine for the back-symmetrisation .
    ! the gaunt coeffients are symmetric too (since the are calculated
    ! using the real spherical harmonics) . that is why the lm2- and
    ! the lm02- loops are only only going up to lm1 or lm01 and the
    ! summands are multiplied by a factor of 2 in the case of lm1 .ne.
    ! lm2 or lm01 .ne. lm02 .

    ! (see notes by b.drittler)

    ! b.drittler   aug. 1988
    ! -----------------------------------------------------------------------
    use :: mod_datatypes, only: dp
    use :: global_variables
    implicit none
    ! lm-dos
    ! ..
    complex (kind=dp) :: df, ek
    integer :: iend, imt1, ipan1, nsra, irmin, irmax
    ! .. Local Scalars ..
    ! ..
    complex (kind=dp) :: cden(irmd, 0:*), cdenns(*), gmat(lmmaxd, lmmaxd), pns(lmmaxd, lmmaxd, irmind:irmd, 2), qnsi(lmmaxd, lmmaxd), qns(lmmaxd, lmmaxd, irmind:irmd, 2), &
      cdenlm(irmd, *), cwr(irmd, lmmaxd, lmmaxd) ! .. Local Arrays ..
    real (kind=dp) :: cleb(*), rho2ns(irmd, lmpotd), thetas(irid, nfund)
    integer :: icleb(ncleb, 4), ifunm(*), lmsp(*)
    ! ..
    complex (kind=dp) :: cltdf, cone, czero
    real (kind=dp) :: c0ll
    integer :: i, ifun, ir, j, l1, lm1, lm2, lm3, m1
    ! ..
    complex (kind=dp) :: wr(lmmaxd, lmmaxd, irmind:irmd), wr2(lmmaxd, lmmaxd, irmind:irmd)
    ! ..

    data czero/(0.0e0_dp, 0.0e0_dp)/
    data cone/(1.0e0_dp, 0.0e0_dp)/
    logical :: opt

    ! ---> initialize array for complex charge density

    c0ll = 1.0e0_dp/sqrt(16.0e0_dp*atan(1.0e0_dp))
    ! ------------------------------------------------------------------

    ! ---> set up array ek*qns(lm1,lm2) + { gmat(lm3,lm2)*pns(lm1,lm3) }
    ! summed over lm3
    cden(1:irmd, 0:lmaxd) = czero
    cwr(:, :, :) = czero
    ! ---> set up of wr(lm1,lm2) = { pns(lm1,lm3)*qns(lm2,lm3) }
    ! summed over lm3

    ! LM2
    ! LM2
    ! LM1
    do ir = irmin + 1, irmax
      do lm1 = 1, lmmaxd
        do lm2 = 1, lmmaxd
          qnsi(lm1, lm2) = qns(lm1, lm2, ir, 1)
        end do
      end do
      call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, cone, pns(1,1,ir,1), lmmaxd, gmat, lmmaxd, ek, qnsi, lmmaxd)
      call zgemm('N', 'T', lmmaxd, lmmaxd, lmmaxd, cone, pns(1,1,ir,1), lmmaxd, qnsi, lmmaxd, czero, wr(1,1,ir), lmmaxd)
      if (nsra==2) then
        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            qnsi(lm1, lm2) = qns(lm1, lm2, ir, 2)
          end do
        end do
        call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, cone, pns(1,1,ir,2), lmmaxd, gmat, lmmaxd, ek, qnsi, lmmaxd)
        call zgemm('N', 'T', lmmaxd, lmmaxd, lmmaxd, cone, pns(1,1,ir,2), lmmaxd, qnsi, lmmaxd, cone, wr(1,1,ir), lmmaxd)
      end if
      ! IR
      do lm1 = 1, lmmaxd
        do lm2 = 1, lm1 - 1
          wr(lm1, lm2, ir) = wr(lm1, lm2, ir) + wr(lm2, lm1, ir)
        end do
        do lm2 = 1, lmmaxd
          wr2(lm1, lm2, ir) = wr(lm1, lm2, ir)
        end do                     ! ---> first calculate only the spherically
        ! symmetric contribution
      end do
    end do
    ! ---> fill array for complex density of states

    ! lm-dos
    do l1 = 0, lmaxd
      do m1 = -l1, l1
        lm1 = l1*(l1+1) + m1 + 1
        do ir = irmin + 1, irmax
          ! lmlm-dos
          ! lmlm-dos
          ! lmlm-dos
          cden(ir, l1) = cden(ir, l1) + wr(lm1, lm1, ir)
          cdenlm(ir, lm1) = wr(lm1, lm1, ir) ! IR
          do lm2 = 1, lmmaxd       ! M1
            cwr(ir, lm1, lm2) = wr2(lm1, lm2, ir)
          end do                   ! ---> remember that the gaunt coeffients
          ! for that case are 1/sqrt(4 pi)
        end do
      end do                       ! Implicit integration over energies



      do ir = irmin + 1, irmax
        rho2ns(ir, 1) = rho2ns(ir, 1) + c0ll*aimag(cden(ir,l1)*df)
      end do
      ! lm-dos
      ! lm-dos
      ! lm-dos
      if (ipan1>1) then
        do i = imt1 + 1, irmax
          cden(i, l1) = cden(i, l1)*thetas(i-imt1, 1)*c0ll
          ! lmlm-dos
          do m1 = -l1, l1          ! lmlm-dos
            lm1 = l1*(l1+1) + m1 + 1 ! if LDAU, integrate up to MT
            cdenlm(i, lm1) = cdenlm(i, lm1)*thetas(i-imt1, 1)*c0ll ! LDAU
            do lm2 = 1, lmmaxd     ! LDAU
              cwr(i, lm1, lm2) = cwr(i, lm1, lm2)*thetas(i-imt1, 1)*c0ll ! LDAU
              ! lmlm-dos
              if (opt('LDA+U   ')) then ! lm-dos
                cwr(i, lm1, lm2) = czero
              end if
            end do                 ! L1
          end do

        end do
      end if

    end do                         ! ---> calculate the non spherically
    ! symmetric contribution

    if (ipan1>1) then
      cdenns(1:irmd) = 0.0e0_dp
    end if

    do j = 1, iend
      lm1 = icleb(j, 1)
      lm2 = icleb(j, 2)
      lm3 = icleb(j, 3)
      cltdf = df*cleb(j)
      ! IF (IPAN1.GT.1) THEN


      do ir = irmin + 1, irmax
        rho2ns(ir, lm3) = rho2ns(ir, lm3) + aimag(cltdf*wr(lm1,lm2,ir))
      end do

      if (ipan1>1 .and. lmsp(lm3)>0) then

        ifun = ifunm(lm3)
        do i = imt1 + 1, irmax
          cdenns(i) = cdenns(i) + cleb(j)*wr(lm1, lm2, i)*thetas(i-imt1, ifun)
        end do
        ! Added IRMIN,IRMAX 1.7.2014
      end if
      ! lm-dos
    end do
    ! -----------------------------------------------------------------------

  end subroutine rhoout

end module mod_rhoout
