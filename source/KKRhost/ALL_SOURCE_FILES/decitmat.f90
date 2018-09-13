module mod_decitmat

contains

subroutine decitmat(eryd, zat, ipan, rr, dror, visp, ircut, rirc, krel, nsra, &
  ins, tmatll, loflm, idoldau, lopt, wldauav, solver, soctl, ctl, zrel, vtrel, &
  btrel, drdi, r2drdi, ipand, irmd, lmaxd, lmaxdp1, lm2d, lmmaxd)
  ! **********************************************************************
  ! *                                                                    *
  ! * A modified form of the CALCTMAT routine to deal with the host      *
  ! * t-matrices in case of decimation                                   *
  ! *                                                                    *
  ! * Non-spherical potential not implemented yet, neither LDA+U         *
  ! *                                                                    *
  ! **********************************************************************
   use mod_beshan
  use :: mod_datatypes, only: dp
   use mod_drvreltmat
   use mod_regsol
   use mod_wfmesh
  use mod_cinit
  implicit none

  ! Parameters ..
  real (kind=dp) :: cvlight
  parameter (cvlight=274.0720442e0_dp)
  complex (kind=dp) :: ci
  parameter (ci=(0e0_dp,1e0_dp))

  ! Scalar arguments ..
  integer :: idoldau, ipan, krel, lopt, nsra, ins, zrel
  integer :: ipand, irmd, lm2d, lmaxd, lmaxdp1, lmmaxd
  real (kind=dp) :: zat, rirc, wldauav
  complex (kind=dp) :: eryd
  character (len=10) :: solver

  ! Array arguments ..
  integer :: ircut(0:ipand), loflm(lm2d)
  real (kind=dp) :: rr(irmd), dror(irmd), visp(irmd)
  complex (kind=dp) :: tmatll(lmmaxd, lmmaxd)
  real (kind=dp) :: soctl(krel*lmaxd+1)
  real (kind=dp) :: ctl(krel*lmaxd+1)
  real (kind=dp) :: vtrel(irmd*krel+(1-krel))
  real (kind=dp) :: btrel(irmd*krel+(1-krel))
  real (kind=dp) :: drdi(irmd), r2drdi(irmd*krel+(1-krel))

  ! Local scalars ..
  integer :: ll, lm1
  real (kind=dp) :: rirc1
  complex (kind=dp) :: ek, carg, qf, hlw, blw

  ! Local arrays ..
  real (kind=dp) :: cutoff(irmd)
  real (kind=dp), allocatable :: rs(:, :), s(:)
  complex (kind=dp), allocatable :: bessjw(:), bessyw(:), hankws(:), dlogdp(:)
  complex (kind=dp), allocatable :: tmat(:), mass(:), hamf(:, :), fz(:, :), pz(:, :)



  call cinit(lmmaxd*lmmaxd, tmatll)
  ! ================================================================= KREL
  if (krel==0) then
    allocate (bessjw(0:lmaxdp1), bessyw(0:lmaxdp1), stat=lm1)
    if (lm1/=0) stop '    Allocate BESSJW/BESSYW'
    allocate (hankws(0:lmaxdp1), dlogdp(0:lmaxd), stat=lm1)
    if (lm1/=0) stop '    Allocate HANKWS/DLOGFP'
    allocate (tmat(0:lmaxd), mass(irmd), stat=lm1)
    if (lm1/=0) stop '    Allocate TMAT/MASS'
    allocate (hamf(irmd,0:lmaxd), fz(irmd,0:lmaxd), stat=lm1)
    if (lm1/=0) stop '    Allocate HAMF/FZ'
    allocate (pz(irmd,0:lmaxd), stat=lm1)
    if (lm1/=0) stop '    Allocate PZ'
    allocate (rs(irmd,0:lmaxd), s(0:lmaxd), stat=lm1)
    if (lm1/=0) stop '    Allocate RS/S'
    rirc1 = 1e0_dp/rirc
    call wfmesh(eryd, ek, cvlight, nsra, zat, rr, s, rs, ircut(ipan), irmd, &
      lmaxd)

    carg = rirc*ek
    call beshan(hankws, bessjw, bessyw, carg, lmaxdp1)
    do ll = 0, lmaxdp1
      hankws(ll) = bessyw(ll) - ci*bessjw(ll)
    end do

    call regsol(cvlight, eryd, nsra, dlogdp, fz, hamf, mass, pz, dror, rr, s, &
      visp, zat, ipan, ircut, idoldau, lopt, wldauav, cutoff, irmd, ipand, &
      lmaxd)

    ! ----------------------------------------------------------------------
    ! --> determine KREL=0 t - matrix

    do ll = 0, lmaxd
      qf = real(ll, kind=dp)*rirc1
      hlw = hankws(ll)*dlogdp(ll)
      blw = bessjw(ll)*dlogdp(ll)

      hlw = qf*hankws(ll) - ek*hankws(ll+1) - hlw
      blw = blw - qf*bessjw(ll) + ek*bessjw(ll+1)
      hlw = hlw*ek
      tmat(ll) = blw/hlw
    end do

    ! --> spherical/non-spherical

    if (ins==0) then
      ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      do lm1 = 1, lmmaxd
        tmatll(lm1, lm1) = tmat(loflm(lm1))
      end do
      ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    else
      ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      stop ' not implemented'
      ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    end if
    deallocate (bessjw, bessyw, hankws, dlogdp, stat=lm1)
    if (lm1/=0) stop '    Deallocate'
    deallocate (tmat, mass, hamf, fz, pz, stat=lm1)
    if (lm1/=0) stop '    Deallocate'
    deallocate (rs, s, stat=lm1)
    if (lm1/=0) stop '    Deallocate'
    ! ----------------------------------------------------------------------
  else                             ! KREL
    call drvreltmat(eryd, tmatll, vtrel, btrel, rr, drdi, r2drdi, zrel, &
      ircut(ipan), solver, soctl, ctl, lmmaxd, lmaxd, irmd)
  end if
  ! ================================================================= KREL
end subroutine decitmat

end module mod_decitmat
