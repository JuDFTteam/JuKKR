! Added IRMIN,IRMAX 1.7.2014
    Subroutine rhoout(cden, df, gmat, ek, pns, qns, rho2ns, thetas, ifunm, &
      ipan1, imt1, irmin, irmax, lmsp, cdenns, nsra, cleb, icleb, iend, &
      cdenlm, cwr) ! lm-dos
      Use mod_datatypes, Only: dp
!-----------------------------------------------------------------------
!
!     calculates the charge density from r(irmin) to r(irc)
!      in case of a non spherical input potential .
!
!     fills the array cden for the complex density of states
!
!     attention : the gaunt coeffients are stored in index array
!                   (see subroutine gaunt)
!
!     the structured part of the greens-function (gmat) is symmetric in
!       its lm-indices , therefore only one half of the matrix is
!       calculated in the subroutine for the back-symmetrisation .
!       the gaunt coeffients are symmetric too (since the are calculated
!       using the real spherical harmonics) . that is why the lm2- and
!       the lm02- loops are only only going up to lm1 or lm01 and the
!       summands are multiplied by a factor of 2 in the case of lm1 .ne.
!       lm2 or lm01 .ne. lm02 .
!
!             (see notes by b.drittler)
!
!                               b.drittler   aug. 1988
!-----------------------------------------------------------------------
!     .. Parameters ..
      Implicit None
      Include 'inc.p'
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *********************************************************************
!
!     ..
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
      Integer :: lmmaxd
      Parameter (lmmaxd=(krel+1)*(lmaxd+1)**2)
      Integer :: lmpotd
      Parameter (lmpotd=(lpotd+1)**2)
      Integer :: irmind
      Parameter (irmind=irmd-irnsd)
! lm-dos
!     ..
      Complex (Kind=dp) :: df, ek
      Integer :: iend, imt1, ipan1, nsra, irmin, irmax
!     .. Local Scalars ..
!     ..
      Complex (Kind=dp) :: cden(irmd, 0:*), cdenns(*), gmat(lmmaxd, lmmaxd), &
        pns(lmmaxd, lmmaxd, irmind:irmd, 2), qnsi(lmmaxd, lmmaxd), &
        qns(lmmaxd, lmmaxd, irmind:irmd, 2), cdenlm(irmd, *), &
        cwr(irmd, lmmaxd, lmmaxd) !     .. Local Arrays ..
      Real (Kind=dp) :: cleb(*), rho2ns(irmd, lmpotd), thetas(irid, nfund)
      Integer :: icleb(ncleb, 4), ifunm(*), lmsp(*)
!     ..
!     .. External Subroutines ..
      Complex (Kind=dp) :: cltdf, cone, czero
      Real (Kind=dp) :: c0ll
      Integer :: i, ifun, ir, j, l1, lm1, lm2, lm3, m1
!     ..
!     .. Intrinsic Functions ..
      Complex (Kind=dp) :: wr(lmmaxd, lmmaxd, irmind:irmd), &
        wr2(lmmaxd, lmmaxd, irmind:irmd)
!     ..
!     .. Save statement ..
      External :: zgemm
!     ..
!     .. Data statements ..
      Intrinsic :: atan, aimag, sqrt
!     ..
!
      Save
!     C0LL = 1/sqrt(4*pi)
!
      Data czero/(0.0E0_dp, 0.0E0_dp)/
      Data cone/(1.0E0_dp, 0.0E0_dp)/
      Logical :: opt
!
!---> initialize array for complex charge density
!
      c0ll = 1.0E0_dp/sqrt(16.0E0_dp*atan(1.0E0_dp))
!------------------------------------------------------------------
!
!---> set up array ek*qns(lm1,lm2) + { gmat(lm3,lm2)*pns(lm1,lm3) }
!                                      summed over lm3
      cden(1:irmd, 0:lmaxd) = czero
      cwr(:, :, :) = czero
!---> set up of wr(lm1,lm2) = { pns(lm1,lm3)*qns(lm2,lm3) }
!                                               summed over lm3

! LM2
! LM2
! LM1
      Do ir = irmin + 1, irmax
        Do lm1 = 1, lmmaxd
          Do lm2 = 1, lmmaxd
            qnsi(lm1, lm2) = qns(lm1, lm2, ir, 1)
          End Do
        End Do
        Call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, cone, pns(1,1,ir,1), &
          lmmaxd, gmat, lmmaxd, ek, qnsi, lmmaxd)
        Call zgemm('N', 'T', lmmaxd, lmmaxd, lmmaxd, cone, pns(1,1,ir,1), &
          lmmaxd, qnsi, lmmaxd, czero, wr(1,1,ir), lmmaxd)
        If (nsra==2) Then
          Do lm1 = 1, lmmaxd
            Do lm2 = 1, lmmaxd
              qnsi(lm1, lm2) = qns(lm1, lm2, ir, 2)
            End Do
          End Do
          Call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, cone, pns(1,1,ir,2), &
            lmmaxd, gmat, lmmaxd, ek, qnsi, lmmaxd)
          Call zgemm('N', 'T', lmmaxd, lmmaxd, lmmaxd, cone, pns(1,1,ir,2), &
            lmmaxd, qnsi, lmmaxd, cone, wr(1,1,ir), lmmaxd)
        End If
! IR
        Do lm1 = 1, lmmaxd
          Do lm2 = 1, lm1 - 1
            wr(lm1, lm2, ir) = wr(lm1, lm2, ir) + wr(lm2, lm1, ir)
          End Do !
          Do lm2 = 1, lmmaxd
            wr2(lm1, lm2, ir) = wr(lm1, lm2, ir)
          End Do !---> first calculate only the spherically symmetric contribution
        End Do !
      End Do !
!---> fill array for complex density of states
!
! lm-dos
      Do l1 = 0, lmaxd
        Do m1 = -l1, l1
          lm1 = l1*(l1+1) + m1 + 1
          Do ir = irmin + 1, irmax
! lmlm-dos
! lmlm-dos
! lmlm-dos
            cden(ir, l1) = cden(ir, l1) + wr(lm1, lm1, ir)
            cdenlm(ir, lm1) = wr(lm1, lm1, ir) ! IR
            Do lm2 = 1, lmmaxd ! M1
              cwr(ir, lm1, lm2) = wr2(lm1, lm2, ir) !
            End Do !---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
          End Do !
        End Do ! Implicit integration over energies
!
!
!
        Do ir = irmin + 1, irmax
          rho2ns(ir, 1) = rho2ns(ir, 1) + c0ll*aimag(cden(ir,l1)*df)
        End Do
! lm-dos
! lm-dos
! lm-dos
        If (ipan1>1) Then
          Do i = imt1 + 1, irmax
            cden(i, l1) = cden(i, l1)*thetas(i-imt1, 1)*c0ll
! lmlm-dos
            Do m1 = -l1, l1 ! lmlm-dos
              lm1 = l1*(l1+1) + m1 + 1 ! if LDAU, integrate up to MT
              cdenlm(i, lm1) = cdenlm(i, lm1)*thetas(i-imt1, 1)*c0ll ! LDAU
              Do lm2 = 1, lmmaxd ! LDAU
                cwr(i, lm1, lm2) = cwr(i, lm1, lm2)*thetas(i-imt1, 1)*c0ll ! LDAU
! lmlm-dos
                If (opt('LDA+U   ')) Then ! lm-dos
                  cwr(i, lm1, lm2) = czero
                End If
              End Do ! L1
            End Do !

          End Do
        End If
!
      End Do !---> calculate the non spherically symmetric contribution
!
      If (ipan1>1) Then
        cdenns(1:irmd) = 0.0E0_dp
      End If
!
      Do j = 1, iend
        lm1 = icleb(j, 1)
        lm2 = icleb(j, 2)
        lm3 = icleb(j, 3)
        cltdf = df*cleb(j)
!       IF (IPAN1.GT.1) THEN


        Do ir = irmin + 1, irmax
          rho2ns(ir, lm3) = rho2ns(ir, lm3) + aimag(cltdf*wr(lm1,lm2,ir))
        End Do

        If (ipan1>1 .And. lmsp(lm3)>0) Then

          ifun = ifunm(lm3)
          Do i = imt1 + 1, irmax
            cdenns(i) = cdenns(i) + cleb(j)*wr(lm1, lm2, i)*thetas(i-imt1, &
              ifun)
          End Do
! Added IRMIN,IRMAX 1.7.2014
        End If
! lm-dos
      End Do
!-----------------------------------------------------------------------
!
    End Subroutine
