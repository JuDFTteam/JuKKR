! Added IRMIN 1.7.2014  &
    Subroutine pnstmat(drdi, ek, icst, pz, qz, fz, sz, pns, tmatll, vins, &
      irmin, ipan, ircut, nsra, cleb, icleb, iend, loflm, tmat, lkonv, &
      idoldau, lopt, lmlo, lmhi, wldau, wldauav, cutoff, alpha0) ! LLY
      Use mod_datatypes, Only: dp
      Implicit None
!     .. Parameters ..
      Include 'inc.p'
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *  LDA+U implementation     Mar. 2002-Dec.2004                      *
! *                           ph.mavropoulos, h. ebert, v. popescu    *
! *                                                                   *
! *********************************************************************
!..
!.. Scalar Arguments ..
!..
!.. Array Arguments ..
      Integer :: irmind
      Parameter (irmind=irmd-irnsd)
      Integer :: mmaxd
      Parameter (mmaxd=2*lmaxd+1)
      Integer :: lmmaxd
      Parameter (lmmaxd=(krel+1)*(lmaxd+1)**2)
      Integer :: lmpotd
      Parameter (lmpotd=(lpotd+1)**2)
      Complex (Kind=dp) :: czero
      Parameter (czero=(0.E0_dp,0.E0_dp))
! LLY
!..
      Complex (Kind=dp) :: ek
      Integer :: icst, idoldau, iend, ipan, lkonv, lopt, nsra, lmlo, lmhi, &
        irmin
      Real (Kind=dp) :: wldauav
!.. Local Scalars ..
!..
      Complex (Kind=dp) :: fz(irmd, 0:lmaxd), pns(lmmaxd, lmmaxd, irmind:irmd, &
        2), pz(irmd, 0:lmaxd), qz(irmd, 0:lmaxd), sz(irmd, 0:lmaxd), &
        tmat(0:lmaxd), tmatll(lmmaxd, lmmaxd), alpha0(lmmaxd, lmmaxd) !.. Local Arrays ..
      Real (Kind=dp) :: cleb(ncleb, 2), drdi(irmd), vins(irmind:irmd, lmpotd)
      Real (Kind=dp) :: wldau(mmaxd, mmaxd), cutoff(irmd)
      Integer :: icleb(ncleb, 4), ircut(0:ipand), loflm(*)
!..
!.. External Subroutines ..
      Integer :: i, ir, lm1, lm2, lmmkonv, m1, m2, irmax


      Complex (Kind=dp) :: ar(lmmaxd, lmmaxd), cmat(lmmaxd, lmmaxd, irmind: &
        irmd), dmat(lmmaxd, lmmaxd, irmind:irmd), efac(lmmaxd), &
        pzekdr(lmmaxd, irmind:irmd, 2), pzlm(lmmaxd, irmind:irmd, 2), &
        qzekdr(lmmaxd, irmind:irmd, 2), qzlm(lmmaxd, irmind:irmd, 2)
      Real (Kind=dp) :: vnspll(lmmaxd, lmmaxd, irmind:irmd)
! ======================================================================
! LDA+U
      External :: regns, vllns, wftsca, zgemm
! Add WLDAU to non-spherical porential VINS in case of LDA+U
      irmax = ircut(ipan)
! Use the average wldau (=wldauav) and calculate the deviation
      Call vllns(vnspll, vins, cleb, icleb, iend, irmd, ncleb, lmpotd, irmind, &
        lmmaxd)
      If (lkonv/=lmaxd) Then
        lmmkonv = (lkonv+1)*(lkonv+1)
        Do lm1 = 1, lmmaxd
          Do lm2 = lmmkonv + 1, lmmaxd
            Do i = irmind, irmd
              vnspll(lm2, lm1, i) = 0.0E0_dp
              vnspll(lm1, lm2, i) = 0.0E0_dp
            End Do
          End Do
        End Do
      Else
        lmmkonv = lmmaxd
      End If
! of wldau from this. Use the deviation in the Born series
! for the non-spherical wavefunction, while the average is
! used for the spherical wavefunction.


! -> First add wldau to all elements ...


      If (idoldau==1 .And. lopt>=0) Then
        Do ir = irmind, irmd
! ... and then subtract average from diag. elements


          Do lm2 = lmlo, lmhi
            m2 = lm2 - lmlo + 1
            Do lm1 = lmlo, lmhi
              m1 = lm1 - lmlo + 1
              vnspll(lm1, lm2, ir) = vnspll(lm1, lm2, ir) + &
                wldau(m1, m2)*cutoff(ir)
            End Do
! LDA+U
! ======================================================================

            vnspll(lm2, lm2, ir) = vnspll(lm2, lm2, ir) - wldauav*cutoff(ir)
          End Do
        End Do
      End If
!---> get wfts of same magnitude by scaling with efac

! Added IRMIN,IRMAX 1.7.2014
      pzlm(:, irmind:irmd, :) = czero
      qzlm(:, irmind:irmd, :) = czero
      pzekdr(:, irmind:irmd, :) = czero
      qzekdr(:, irmind:irmd, :) = czero
      cmat(:, :, irmind:irmd) = czero
      dmat(:, :, irmind:irmd) = czero

!---> determine the regular non sph. wavefunction

      Call wftsca(drdi, efac, pz, qz, fz, sz, nsra, pzlm, qzlm, pzekdr, &
        qzekdr, ek, loflm, irmind, irmd, irmin, irmax, lmaxd, lmmaxd) ! Added IRMIN,IRMAX 1.7.2014  &



      Call regns(ar, tmatll, efac, pns, vnspll, icst, ipan, ircut, pzlm, qzlm, &
        pzekdr, qzekdr, ek, pns(1,1,irmind,1), cmat, pns(1,1,irmind,2), dmat, &
        nsra, irmind, irmd, irmin, irmax, ipand, lmmaxd)
! LLY non-spher. contribution to alpha matrix
! LLY Drittler PhD eq. 3.106
      Do lm1 = 1, lmmkonv
        tmatll(lm1, lm1) = tmatll(lm1, lm1) + tmat(loflm(lm1))
      End Do
! on output.

      Do lm2 = 1, lmmkonv
        Do lm1 = 1, lmmkonv
          ar(lm1, lm2) = alpha0(lm1, lm1)*ar(lm1, lm2)
        End Do ! Added IRMIN 1.7.2014  &
      End Do
      alpha0(1:lmmaxd, 1:lmmaxd) = ar(1:lmmaxd, 1:lmmaxd) ! LLY
!     .. Parameters ..
! set to 1 if NEWSOSOL under RUNOPT, otherwise 0
      Return
    End Subroutine
