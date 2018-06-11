    Subroutine phicalc(iatom, lphi, visp, ipan, ircut, r, drdi, z, erefldau, &
      phi, nspin, nsra)
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

      Use mod_datatypes
      Implicit None
      Include 'inc.p'
      Integer :: npotd
      Parameter (npotd=(2*krel+(1-krel)*nspind)*natypd)
      Real (Kind=dp) :: cvlight
      Parameter (cvlight=274.0720442E0_dp)

      Complex (Kind=dp) :: phi(irmd), pz(irmd, 0:lmaxd), fz(irmd, 0:lmaxd)
      Complex (Kind=dp) :: hamf(irmd, 0:lmaxd), mass(irmd), dlogdp(0:lmaxd)

      Real (Kind=dp) :: r(irmd, natypd), visp(irmd, npotd), z(natypd)
      Real (Kind=dp) :: drdi(irmd, natypd)
      Real (Kind=dp) :: rs(irmd, 0:lmaxd), s(0:lmaxd), dror(irmd), vpot(irmd)
      Real (Kind=dp) :: erefldau

      Integer :: ipan(natypd), ircut(0:ipand, natypd)
      Integer :: iatom, ipot1, nsra, nspin, irc1, ipan1
! points at the correct potential,
      Integer :: lphi
! i.e., 1,3,5,.. for NSPIN=2.
      Integer :: ir, l1
      Real (Kind=dp) :: wint(irmd), wnorm, cutoff(irmd)
      Complex (Kind=dp) :: cnorm, ez, czero

      czero = cmplx(0.E0_dp, 0.E0_dp, kind=dp)
      ipan1 = ipan(iatom)
      irc1 = ircut(ipan1, iatom)
! -> set VPOT = [ V(UP) + V(DN) ]/2  for IATOM
      ipot1 = (iatom-1)*nspin + 1 !    only for the spherical part.


! -> Prepare and call REGSOL


      Call dcopy(irmd, visp(1,ipot1), 1, vpot(1), 1)
      If (nspin==2) Then
        Call daxpy(irmd, 1.E0_dp, visp(1,ipot1+1), 1, vpot(1), 1)
        Call dscal(irmd, 0.5E0_dp, vpot, 1)
      End If

! --> this call of regsol within a non-lda+u calculation

      Do l1 = 0, lmaxd
        If (nsra==2) Then
          s(l1) = sqrt(real(l1*l1+l1+1,kind=dp)-4.0E0_dp*z(iatom)*z(iatom)/( &
            cvlight*cvlight))
          If (z(iatom)==0.0E0_dp) s(l1) = real(l1, kind=dp)
        Else
          s(l1) = real(l1, kind=dp)
        End If
        rs(1, l1) = 0.0E0_dp
        Do ir = 2, irmd
          rs(ir, l1) = r(ir, iatom)**s(l1)
        End Do
      End Do

      Do ir = 2, irc1
        dror(ir) = drdi(ir, iatom)/r(ir, iatom)
      End Do
      ez = cmplx(erefldau, 0.E0_dp, kind=dp)
! The result of regsol is (R*r)/r**l (in non-sra and similar in sra)


      Call regsol(cvlight, ez, nsra, dlogdp, fz, hamf, mass, pz, dror, &
        r(1,iatom), s, vpot, z(iatom), ipan(iatom), ircut(0,iatom), 0, -1, &
        0E0_dp, cutoff, irmd, ipand, lmaxd)
! --> set current angular momentum as LPHI


! -> Copy result to PHI (only large component)


      If (nsra==2) Then
        Do ir = 2, irc1
          pz(ir, lphi) = pz(ir, lphi)*rs(ir, lphi)
          fz(ir, lphi) = fz(ir, lphi)/cvlight*rs(ir, lphi)/cvlight
        End Do
      Else
        Do ir = 2, irc1
          pz(ir, lphi) = pz(ir, lphi)*rs(ir, lphi)
          fz(ir, lphi) = czero
        End Do
      End If
! -> Normalise in sphere:
!    Note: If we have normalised in cell instead of sphere, the choice
!    M=0 would be arbitrary. One would have different normalisation for
      Call zcopy(irmd, pz(1,lphi), 1, phi(1), 1)
!    each M, if the cell has no cubic symmetry.  Here we normalise in
!    sphere, to avoid this inconsistency.



! --> integrate, normalisation factor = WNORM

      Do ir = 1, irc1
        wint(ir) = real(conjg(phi(ir))*phi(ir))
      End Do

! --> normalise PZ,FZ to unit probability in WS cell


      Call simpk(wint, wnorm, ipan1, ircut(0,iatom), drdi(1,iatom))
! *********************************************************************
! *                                                                   *
! * Calculates test functions PHI for LDA+U.                          *
      cnorm = 1.E0_dp/sqrt(wnorm)
      Call zscal(irc1, cnorm, phi, 1)
! * Only spherical part of the potential is used.                     *
      Return
    End Subroutine
