!------------------------------------------------------------------------------------
!> Summary: Auxiliary function to calculate the single site t-matrix for LDA+U implementation 
!> Author: Ph. Mavropoulos, H. Ebert, V. Popescu
!> Auxiliary function to calculate the single site t-matrix for LDA+U implementation 
!------------------------------------------------------------------------------------
module mod_pnstmat

contains

  !-------------------------------------------------------------------------------
  !> Summary: Auxiliary function to calculate the single site t-matrix for LDA+U implementation 
  !> Author: Ph. Mavropoulos, H. Ebert, V. Popescu
  !> Category: lda+u, single-site, KKRhost 
  !> Deprecated: False 
  !> Auxiliary function to calculate the single site t-matrix for LDA+U implementation 
  !-------------------------------------------------------------------------------
  subroutine pnstmat(drdi,ek,icst,pz,qz,fz,sz,pns,tmatll,vins,irmin,ipan,ircut,nsra,& 
    cleb,icleb,iend,loflm,tmat,lkonv,idoldau,lopt,lmlo,lmhi,wldau,wldauav,cutoff,   &
    alpha0)                ! LLY

    use :: mod_datatypes, only: dp
    use :: global_variables
    use :: mod_vllns
    use :: mod_wftsca
    use :: mod_regns
    use :: mod_constants, only: czero
    implicit none
    ! LLY
    ! ..
    complex (kind=dp) :: ek
    integer :: icst, idoldau, iend, ipan, lkonv, lopt, nsra, lmlo, lmhi, irmin
    real (kind=dp) :: wldauav
    ! .. Local Scalars ..
    ! ..
    complex (kind=dp) :: fz(irmd, 0:lmaxd), pns(lmmaxd, lmmaxd, irmind:irmd, 2), pz(irmd, 0:lmaxd), qz(irmd, 0:lmaxd), sz(irmd, 0:lmaxd), tmat(0:lmaxd), tmatll(lmmaxd, lmmaxd), &
      alpha0(lmmaxd, lmmaxd)       ! .. Local Arrays ..
    real (kind=dp) :: cleb(ncleb, 2), drdi(irmd), vins(irmind:irmd, lmpotd)
    real (kind=dp) :: wldau(mmaxd, mmaxd), cutoff(irmd)
    integer :: icleb(ncleb, 4), ircut(0:ipand), loflm(*)
    ! ..
    integer :: i, ir, lm1, lm2, lmmkonv, m1, m2, irmax


    complex (kind=dp) :: ar(lmmaxd, lmmaxd), cmat(lmmaxd, lmmaxd, irmind:irmd), dmat(lmmaxd, lmmaxd, irmind:irmd), efac(lmmaxd), pzekdr(lmmaxd, irmind:irmd, 2), &
      pzlm(lmmaxd, irmind:irmd, 2), qzekdr(lmmaxd, irmind:irmd, 2), qzlm(lmmaxd, irmind:irmd, 2)
    real (kind=dp) :: vnspll(lmmaxd, lmmaxd, irmind:irmd)
    ! ======================================================================
    ! LDA+U
    ! Add WLDAU to non-spherical porential VINS in case of LDA+U
    irmax = ircut(ipan)
    ! Use the average wldau (=wldauav) and calculate the deviation
    call vllns(vnspll, vins, cleb, icleb, iend, irmd, ncleb, lmpotd, irmind, lmmaxd)
    if (lkonv/=lmaxd) then
      lmmkonv = (lkonv+1)*(lkonv+1)
      do lm1 = 1, lmmaxd
        do lm2 = lmmkonv + 1, lmmaxd
          do i = irmind, irmd
            vnspll(lm2, lm1, i) = 0.0e0_dp
            vnspll(lm1, lm2, i) = 0.0e0_dp
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
            vnspll(lm1, lm2, ir) = vnspll(lm1, lm2, ir) + wldau(m1, m2)*cutoff(ir)
          end do
          ! LDA+U
          ! ======================================================================

          vnspll(lm2, lm2, ir) = vnspll(lm2, lm2, ir) - wldauav*cutoff(ir)
        end do
      end do
    end if
    ! ---> get wfts of same magnitude by scaling with efac

    ! Added IRMIN,IRMAX 1.7.2014
    pzlm(:, irmind:irmd, :) = czero
    qzlm(:, irmind:irmd, :) = czero
    pzekdr(:, irmind:irmd, :) = czero
    qzekdr(:, irmind:irmd, :) = czero
    cmat(:, :, irmind:irmd) = czero
    dmat(:, :, irmind:irmd) = czero

    ! ---> determine the regular non sph. wavefunction

    call wftsca(drdi,efac,pz,qz,fz,sz,nsra,pzlm,qzlm,pzekdr,qzekdr,ek,loflm,irmind, &
      irmd,irmin,irmax,lmaxd,lmmaxd) ! Added IRMIN,IRMAX
    ! 1.7.2014  &

    call regns(ar,tmatll,efac,pns,vnspll,icst,ipan,ircut,pzlm,qzlm,pzekdr,qzekdr,ek,& 
      pns(1,1,irmind,1),cmat,pns(1,1,irmind,2),dmat,nsra,irmind,irmd,irmin,irmax,   &
      ipand,lmmaxd)
    ! LLY non-spher. contribution to alpha matrix
    ! LLY Drittler PhD eq. 3.106
    do lm1 = 1, lmmkonv
      tmatll(lm1, lm1) = tmatll(lm1, lm1) + tmat(loflm(lm1))
    end do
    ! on output.

    do lm2 = 1, lmmkonv
      do lm1 = 1, lmmkonv
        ar(lm1, lm2) = alpha0(lm1, lm1)*ar(lm1, lm2)
      end do                       ! Added IRMIN 1.7.2014  &
    end do
    alpha0(1:lmmaxd, 1:lmmaxd) = ar(1:lmmaxd, 1:lmmaxd) ! LLY
    ! .. Parameters ..
    ! set to 1 if NEWSOSOL under RUNOPT, otherwise 0
    return
  end subroutine pnstmat

end module mod_pnstmat
