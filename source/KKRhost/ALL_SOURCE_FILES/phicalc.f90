!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Calculates test functions PHI for LDA+U
!> Author: Ph. Mavropoulos
!> Calculates test functions PHI for LDA+U. Only spherical part of the potential is used.
!> PHI are then normalised within Voronoi Circumscribing Sphere.
!> In SRA treatment, only large component is considered as important and normalised 
!> although small component is also calculated. `PHI` is here one-dimensional array 
!> -- `PHI(IRMD)` -- so the subroutine must be called for each atom seperately.
!------------------------------------------------------------------------------------
module mod_phicalc

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculates test functions PHI for LDA+U
  !> Author: Ph. Mavropoulos
  !> Category: lda+u, KKRhost 
  !> Deprecated: False 
  !> Calculates test functions PHI for LDA+U. Only spherical part of the potential is used.
  !> PHI are then normalised within Voronoi Circumscribing Sphere.
  !> In SRA treatment, only large component is considered as important and normalised 
  !> although small component is also calculated. `PHI` is here one-dimensional array 
  !> -- `PHI(IRMD)` -- so the subroutine must be called for each atom seperately.
  !-------------------------------------------------------------------------------
  subroutine phicalc(iatom,lphi,visp,ipan,ircut,r,drdi,z,erefldau,phi,nspin,nsra)

    use :: mod_datatypes
    use :: global_variables
    use :: mod_regsol
    use :: mod_simpk
    use :: mod_constants, only: cvlight, czero
    implicit none
    real (kind=dp), parameter :: eps = 1.0d-12
    ! .. Input variables
    integer, intent(in) :: lphi  !! points at the correct potential, i.e., 1,3,5,.. for NSPIN=2.
    integer, intent(in) :: nsra
    integer, intent(in) :: nspin !! Counter for spin directions
    integer, intent(in) :: iatom
    real (kind=dp), intent(in) :: erefldau !! the energies of the projector's wave functions (REAL)
    integer, dimension(natypd), intent(in) :: ipan  !! Number of panels in non-MT-region
    integer, dimension(0:ipand, natypd), intent(in) :: ircut  !! R points of panel borders
    real (kind=dp), dimension(natypd), intent(in) :: z !! Nuclear charge
    real (kind=dp), dimension(irmd, natypd), intent(in) :: r    !! Set of real space vectors (in a.u.)
    real (kind=dp), dimension(irmd, npotd), intent(in)  :: visp !! Spherical part of the potential
    real (kind=dp), dimension(irmd, natypd), intent(in) :: drdi !! Derivative dr/di
    ! .. In/Out variables
    complex (kind=dp), dimension(irmd), intent(inout) :: phi !! Test function
    ! .. Local variables
    integer :: ir, l1, ipot1, irc1, ipan1
    real (kind=dp) :: wnorm
    complex (kind=dp) :: cnorm, ez
    real (kind=dp), dimension(0:lmaxd)  :: s
    real (kind=dp), dimension(irmd)     :: dror
    real (kind=dp), dimension(irmd)     :: vpot
    real (kind=dp), dimension(irmd)     :: wint
    real (kind=dp), dimension(irmd)     :: cutoff
    real (kind=dp), dimension(irmd, 0:lmaxd) :: rs
    complex (kind=dp), dimension(irmd)    :: mass
    complex (kind=dp), dimension(0:lmaxd) :: dlogdp
    complex (kind=dp), dimension(irmd, 0:lmaxd) :: pz
    complex (kind=dp), dimension(irmd, 0:lmaxd) :: fz
    complex (kind=dp), dimension(irmd, 0:lmaxd) :: hamf

    ipan1 = ipan(iatom)
    irc1 = ircut(ipan1, iatom)
    ! -> set VPOT = [ V(UP) + V(DN) ]/2  for IATOM
    ipot1 = (iatom-1)*nspin + 1    ! only for the spherical part.

    ! -> Prepare and call REGSOL
    call dcopy(irmd, visp(1,ipot1), 1, vpot(1), 1)
    if (nspin==2) then
      call daxpy(irmd, 1.e0_dp, visp(1,ipot1+1), 1, vpot(1), 1)
      call dscal(irmd, 0.5e0_dp, vpot, 1)
    end if

    ! --> this call of regsol within a non-lda+u calculation
    do l1 = 0, lmaxd
      if (nsra==2) then
        s(l1) = sqrt(real(l1*l1+l1+1,kind=dp)-4.0e0_dp*z(iatom)*z(iatom)/(cvlight*cvlight))
        if (abs(z(iatom))<eps) s(l1) = real(l1, kind=dp)
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

    call regsol(cvlight,ez,nsra,dlogdp,fz,hamf,mass,pz,dror,r(1,iatom),s,vpot,      &
      z(iatom),ipan(iatom),ircut(0,iatom),0,-1,0e0_dp,cutoff,irmd,ipand,lmaxd)
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

end module mod_phicalc
