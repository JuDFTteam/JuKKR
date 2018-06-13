subroutine phicalc(iatom, lphi, visp, ipan, ircut, r, drdi, z, erefldau, phi, &
  nspin, nsra)
  ! *********************************************************************
  ! *                                                                   *
  ! * Calculates test functions PHI for LDA+U.                          *
  ! * Only spherical part of the potential is used.                     *
  ! * PHI are then normalised within Voronoi Circumscribing Sphere.     *
  ! * In SRA treatment, only large component is considered as important *
  ! * and normalised although small component is also calculated.       *
  ! *                                                                   *
  ! * PHI is here one-dime array -- PHI(IRMD) -- so the subroutine must *
  ! * be called for each atom seperately.                               *
  ! *                                  ph. mavropoulos, juelich 2004    *
  ! *                                                                   *
  ! *********************************************************************

  use :: mod_datatypes
  use global_variables
  implicit none
  real (kind=dp) :: cvlight
  parameter (cvlight=274.0720442e0_dp)

  complex (kind=dp) :: phi(irmd), pz(irmd, 0:lmaxd), fz(irmd, 0:lmaxd)
  complex (kind=dp) :: hamf(irmd, 0:lmaxd), mass(irmd), dlogdp(0:lmaxd)

  real (kind=dp) :: r(irmd, natypd), visp(irmd, npotd), z(natypd)
  real (kind=dp) :: drdi(irmd, natypd)
  real (kind=dp) :: rs(irmd, 0:lmaxd), s(0:lmaxd), dror(irmd), vpot(irmd)
  real (kind=dp) :: erefldau

  integer :: ipan(natypd), ircut(0:ipand, natypd)
  integer :: iatom, ipot1, nsra, nspin, irc1, ipan1
  ! points at the correct potential,
  integer :: lphi
  ! i.e., 1,3,5,.. for NSPIN=2.
  integer :: ir, l1
  real (kind=dp) :: wint(irmd), wnorm, cutoff(irmd)
  complex (kind=dp) :: cnorm, ez, czero

  czero = cmplx(0.e0_dp, 0.e0_dp, kind=dp)
  ipan1 = ipan(iatom)
  irc1 = ircut(ipan1, iatom)
  ! -> set VPOT = [ V(UP) + V(DN) ]/2  for IATOM
  ipot1 = (iatom-1)*nspin + 1      ! only for the spherical part.


  ! -> Prepare and call REGSOL


  call dcopy(irmd, visp(1,ipot1), 1, vpot(1), 1)
  if (nspin==2) then
    call daxpy(irmd, 1.e0_dp, visp(1,ipot1+1), 1, vpot(1), 1)
    call dscal(irmd, 0.5e0_dp, vpot, 1)
  end if

  ! --> this call of regsol within a non-lda+u calculation

  do l1 = 0, lmaxd
    if (nsra==2) then
      s(l1) = sqrt(real(l1*l1+l1+1,kind=dp)-4.0e0_dp*z(iatom)*z(iatom)/( &
        cvlight*cvlight))
      if (z(iatom)==0.0e0_dp) s(l1) = real(l1, kind=dp)
    else
      s(l1) = real(l1, kind=dp)
    end if
    rs(1, l1) = 0.0e0_dp
    do ir = 2, irmd
      rs(ir, l1) = r(ir, iatom)**s(l1)
    end do
  end do

  do ir = 2, irc1
    dror(ir) = drdi(ir, iatom)/r(ir, iatom)
  end do
  ez = cmplx(erefldau, 0.e0_dp, kind=dp)
  ! The result of regsol is (R*r)/r**l (in non-sra and similar in sra)


  call regsol(cvlight, ez, nsra, dlogdp, fz, hamf, mass, pz, dror, r(1,iatom), &
    s, vpot, z(iatom), ipan(iatom), ircut(0,iatom), 0, -1, 0e0_dp, cutoff, &
    irmd, ipand, lmaxd)
  ! --> set current angular momentum as LPHI


  ! -> Copy result to PHI (only large component)


  if (nsra==2) then
    do ir = 2, irc1
      pz(ir, lphi) = pz(ir, lphi)*rs(ir, lphi)
      fz(ir, lphi) = fz(ir, lphi)/cvlight*rs(ir, lphi)/cvlight
    end do
  else
    do ir = 2, irc1
      pz(ir, lphi) = pz(ir, lphi)*rs(ir, lphi)
      fz(ir, lphi) = czero
    end do
  end if
  ! -> Normalise in sphere:
  ! Note: If we have normalised in cell instead of sphere, the choice
  ! M=0 would be arbitrary. One would have different normalisation for
  call zcopy(irmd, pz(1,lphi), 1, phi(1), 1)
  ! each M, if the cell has no cubic symmetry.  Here we normalise in
  ! sphere, to avoid this inconsistency.



  ! --> integrate, normalisation factor = WNORM

  do ir = 1, irc1
    wint(ir) = real(conjg(phi(ir))*phi(ir))
  end do

  ! --> normalise PZ,FZ to unit probability in WS cell


  call simpk(wint, wnorm, ipan1, ircut(0,iatom), drdi(1,iatom))
  ! *********************************************************************
  ! *                                                                   *
  ! * Calculates test functions PHI for LDA+U.                          *
  cnorm = 1.e0_dp/sqrt(wnorm)
  call zscal(irc1, cnorm, phi, 1)
  ! * Only spherical part of the potential is used.                     *
  return
end subroutine phicalc
