! -------------------------------------------------------------------------------
! SUBROUTINE: PNSQNS
! > @note
! > - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to
! Fortran90
! -------------------------------------------------------------------------------
! Added IRMIN 1.7.2014
subroutine pnsqns(ar, cr, dr, drdi, ek, icst, pz, qz, fz, sz, pns, qns, nsra, &
  vins, ipan, irmin, ircut, cleb, icleb, iend, loflm, lkonv, idoldau, lopt, &
  lmlo, lmhi, wldau, wldauav, cutoff, lmax)

  use :: global_variables
  use :: mod_datatypes, only: dp

  implicit none

  ! .. Input variables
  integer, intent (in) :: lmax     ! < Maximum l component in wave function
                                   ! expansion
  integer, intent (in) :: icst     ! < Number of Born approximation
  integer, intent (in) :: iend     ! < Number of nonzero gaunt coefficients
  integer, intent (in) :: ipan     ! < Number of panels in non-MT-region
  integer, intent (in) :: lmlo
  integer, intent (in) :: lmhi
  integer, intent (in) :: lopt     ! < angular momentum QNUM for the atoms on
                                   ! which LDA+U should be applied (-1 to
                                   ! switch it OFF)
  integer, intent (in) :: nsra
  integer, intent (in) :: lkonv
  integer, intent (in) :: irmin    ! < Max R for spherical treatment
  integer, intent (in) :: idoldau  ! < flag to perform LDA+U
  real (kind=dp), intent (in) :: wldauav
  complex (kind=dp), intent (in) :: ek
  integer, dimension (0:ipand), intent (in) :: ircut ! < R points of panel
                                                     ! borders
  integer, dimension (*), intent (in) :: loflm ! < l of lm=(l,m) (GAUNT)
  integer, dimension (ncleb, 4), intent (in) :: icleb ! < Pointer array
  real (kind=dp), dimension (irmd), intent (in) :: drdi ! < Derivative dr/di
  real (kind=dp), dimension (irmd), intent (in) :: cutoff
  real (kind=dp), dimension (ncleb, 2), intent (in) :: cleb ! < GAUNT
                                                            ! coefficients
                                                            ! (GAUNT)
  real (kind=dp), dimension (irmind:irmd, lmpotd), intent (in) :: vins ! <
                                                                       ! Non-spherical
                                                                       ! part
                                                                       ! of
                                                                       ! the
                                                                       ! potential
  real (kind=dp), dimension (mmaxd, mmaxd), intent (in) :: wldau ! < potential
                                                                 ! matrix
  complex (kind=dp), dimension (irmd, 0:lmax), intent (in) :: fz
  complex (kind=dp), dimension (irmd, 0:lmax), intent (in) :: qz
  complex (kind=dp), dimension (irmd, 0:lmax), intent (in) :: sz
  complex (kind=dp), dimension (irmd, 0:lmax), intent (in) :: pz
  complex (kind=dp), dimension (lmmaxd, lmmaxd), intent (in) :: dr
  complex (kind=dp), dimension (lmmaxd, lmmaxd), intent (in) :: ar
  complex (kind=dp), dimension (lmmaxd, lmmaxd), intent (in) :: cr
  complex (kind=dp), dimension (lmmaxd, lmmaxd, irmind:irmd, 2), &
    intent (in) :: pns
  complex (kind=dp), dimension (lmmaxd, lmmaxd, irmind:irmd, 2), &
    intent (in) :: qns
  ! .. Local Scalars
  integer :: i, lm1, lm2, lmmkonv, m1, m2, ir, irmax
  ! .. Local Arrays
  real (kind=dp), dimension (lmmaxd, lmmaxd, irmind:irmd) :: vnspll
  complex (kind=dp), dimension (lmmaxd) :: efac
  complex (kind=dp), dimension (lmmaxd, lmmaxd) :: tmatll
  complex (kind=dp), dimension (lmmaxd, lmmaxd, irmind:irmd) :: dmat
  complex (kind=dp), dimension (lmmaxd, lmmaxd, irmind:irmd) :: cmat
  complex (kind=dp), dimension (lmmaxd, irmind:irmd, 2) :: pzlm
  complex (kind=dp), dimension (lmmaxd, irmind:irmd, 2) :: qzlm
  complex (kind=dp), dimension (lmmaxd, irmind:irmd, 2) :: pzekdr
  complex (kind=dp), dimension (lmmaxd, irmind:irmd, 2) :: qzekdr
  ! .. External Subroutines
  external :: irwns, regns, vllns, wftsca

  irmax = ircut(ipan)              ! Added IRMAX 1.7.2014
  call vllns(vnspll, vins, cleb, icleb, iend, irmd, ncleb, lmpotd, irmind, &
    lmmaxd)
  if (lkonv/=lmax) then
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
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! LDA+U
  ! Add WLDAU to non-spherical porential VINS in case of LDA+U
  ! Use the average wldau (=wldauav) and the deviation
  ! of wldau from this. Use the deviation in the Born series
  ! for the non-spherical wavefunction, while the average is
  ! used for the spherical wavefunction.
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (idoldau==1 .and. lopt>=0) then
    do ir = irmind, irmd
      ! ----------------------------------------------------------------------
      ! First add wldau to all elements
      ! ----------------------------------------------------------------------
      do lm2 = lmlo, lmhi
        m2 = lm2 - lmlo + 1
        do lm1 = lmlo, lmhi
          m1 = lm1 - lmlo + 1
          vnspll(lm1, lm2, ir) = vnspll(lm1, lm2, ir) + &
            wldau(m1, m2)*cutoff(ir)
        end do
        ! -------------------------------------------------------------------
        ! and then subtract average from diag. elements
        ! -------------------------------------------------------------------
        vnspll(lm2, lm2, ir) = vnspll(lm2, lm2, ir) - wldauav*cutoff(ir)
      end do
    end do
  end if
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! LDA+U
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ----------------------------------------------------------------------------
  ! Get wfts of same magnitude by scaling with efac
  ! ----------------------------------------------------------------------------
  call wftsca(drdi, efac, pz, qz, fz, sz, nsra, pzlm, qzlm, pzekdr, qzekdr, &
    ek, loflm, irmind, irmd, irmin, irmax, lmax, lmmaxd) ! Added IRMIN,IRMAX
                                                         ! 1.7.2014
  ! ----------------------------------------------------------------------------
  ! Determine the irregular non sph. wavefunction
  ! ----------------------------------------------------------------------------
  ! Added IRMIN,IRMAX 1.7.2014
  call irwns(cr, dr, efac, qns, vnspll, icst, ipan, ircut, nsra, pzlm, qzlm, &
    pzekdr, qzekdr, qns(1,1,irmind,1), cmat, qns(1,1,irmind,2), dmat, irmind, &
    irmd, irmin, irmax, ipand, lmmaxd)
  ! ----------------------------------------------------------------------------
  ! Determine the regular non sph. wavefunction
  ! ----------------------------------------------------------------------------
  ! Added IRMIN,IRMAX 1.7.2014
  call regns(ar, tmatll, efac, pns, vnspll, icst, ipan, ircut, pzlm, qzlm, &
    pzekdr, qzekdr, ek, pns(1,1,irmind,1), cmat, pns(1,1,irmind,2), dmat, &
    nsra, irmind, irmd, irmin, irmax, ipand, lmmaxd)

  return

end subroutine pnsqns
