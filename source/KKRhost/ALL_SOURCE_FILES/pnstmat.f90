subroutine pnstmat(drdi, ek, icst, pz, qz, fz, sz, pns, tmatll, vins, irmin, &
  ipan, ircut, nsra, cleb, icleb, iend, loflm, tmat, lkonv, & ! Added IRMIN 1.7.2014  &
  idoldau, lopt, lmlo, lmhi, wldau, wldauav, cutoff, alpha0) ! LLY
  implicit none
!     .. Parameters ..
  include 'inc.p'
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
  integer :: irmind
  parameter (irmind=irmd-irnsd)
  integer :: mmaxd
  parameter (mmaxd=2*lmaxd+1)
  integer :: lmmaxd
  parameter (lmmaxd=(krel+1)*(lmaxd+1)**2)
  integer :: lmpotd
  parameter (lmpotd=(lpotd+1)**2)
  double complex :: czero
  parameter (czero=(0.d0,0.d0))
! LLY
!..
  double complex :: ek
  integer :: icst, idoldau, iend, ipan, lkonv, lopt, nsra, lmlo, lmhi, irmin
  double precision :: wldauav
!.. Local Scalars ..
!..
  double complex :: fz(irmd, 0:lmaxd), pns(lmmaxd, lmmaxd, irmind:irmd, 2), &
    pz(irmd, 0:lmaxd), qz(irmd, 0:lmaxd), sz(irmd, 0:lmaxd), tmat(0:lmaxd), &
    tmatll(lmmaxd, lmmaxd), alpha0(lmmaxd, lmmaxd) !.. Local Arrays ..
  double precision :: cleb(ncleb, 2), drdi(irmd), vins(irmind:irmd, lmpotd)
  double precision :: wldau(mmaxd, mmaxd), cutoff(irmd)
  integer :: icleb(ncleb, 4), ircut(0:ipand), loflm(*)
!..
!.. External Subroutines ..
  integer :: i, ir, lm1, lm2, lmmkonv, m1, m2, irmax


  double complex :: ar(lmmaxd, lmmaxd), cmat(lmmaxd, lmmaxd, irmind:irmd), &
    dmat(lmmaxd, lmmaxd, irmind:irmd), efac(lmmaxd), &
    pzekdr(lmmaxd, irmind:irmd, 2), pzlm(lmmaxd, irmind:irmd, 2), &
    qzekdr(lmmaxd, irmind:irmd, 2), qzlm(lmmaxd, irmind:irmd, 2)
  double precision :: vnspll(lmmaxd, lmmaxd, irmind:irmd)
! ======================================================================
! LDA+U
  external :: regns, vllns, wftsca, zgemm
! Add WLDAU to non-spherical porential VINS in case of LDA+U
  irmax = ircut(ipan)
! Use the average wldau (=wldauav) and calculate the deviation
  call vllns(vnspll, vins, cleb, icleb, iend, irmd, ncleb, lmpotd, irmind, &
    lmmaxd)
  if (lkonv/=lmaxd) then
    lmmkonv = (lkonv+1)*(lkonv+1)
    do lm1 = 1, lmmaxd
      do lm2 = lmmkonv + 1, lmmaxd
        do i = irmind, irmd
          vnspll(lm2, lm1, i) = 0.0d0
          vnspll(lm1, lm2, i) = 0.0d0
        end do
      end do
    end do
  else
    lmmkonv = lmmaxd
  end if
! of wldau from this. Use the deviation in the Born series
! for the non-spherical wavefunction, while the average is
! used for the spherical wavefunction.


! -> First add wldau to all elements ...


  if (idoldau==1 .and. lopt>=0) then
    do ir = irmind, irmd
! ... and then subtract average from diag. elements


      do lm2 = lmlo, lmhi
        m2 = lm2 - lmlo + 1
        do lm1 = lmlo, lmhi
          m1 = lm1 - lmlo + 1
          vnspll(lm1, lm2, ir) = vnspll(lm1, lm2, ir) + &
            wldau(m1, m2)*cutoff(ir)
        end do
! LDA+U
! ======================================================================

        vnspll(lm2, lm2, ir) = vnspll(lm2, lm2, ir) - wldauav*cutoff(ir)
      end do
    end do
  end if
!---> get wfts of same magnitude by scaling with efac

! Added IRMIN,IRMAX 1.7.2014
  pzlm(:, irmind:irmd, :) = czero
  qzlm(:, irmind:irmd, :) = czero
  pzekdr(:, irmind:irmd, :) = czero
  qzekdr(:, irmind:irmd, :) = czero
  cmat(:, :, irmind:irmd) = czero
  dmat(:, :, irmind:irmd) = czero

!---> determine the regular non sph. wavefunction

  call wftsca(drdi, efac, pz, qz, fz, sz, nsra, pzlm, qzlm, pzekdr, qzekdr, &
    ek, loflm, irmind, irmd, irmin, irmax, lmaxd, lmmaxd) ! Added IRMIN,IRMAX 1.7.2014  &



  call regns(ar, tmatll, efac, pns, vnspll, icst, ipan, ircut, pzlm, qzlm, &
    pzekdr, qzekdr, ek, pns(1,1,irmind,1), cmat, pns(1,1,irmind,2), dmat, &
    nsra, irmind, irmd, irmin, irmax, & 
    ipand, lmmaxd)
! LLY non-spher. contribution to alpha matrix
! LLY Drittler PhD eq. 3.106
  do lm1 = 1, lmmkonv
    tmatll(lm1, lm1) = tmatll(lm1, lm1) + tmat(loflm(lm1))
  end do
! on output.

  do lm2 = 1, lmmkonv
    do lm1 = 1, lmmkonv
      ar(lm1, lm2) = alpha0(lm1, lm1)*ar(lm1, lm2) 
    end do ! Added IRMIN 1.7.2014  &
  end do
  alpha0(1:lmmaxd, 1:lmmaxd) = ar(1:lmmaxd, 1:lmmaxd) ! LLY
!     .. Parameters ..
! set to 1 if NEWSOSOL under RUNOPT, otherwise 0
  return
end subroutine
