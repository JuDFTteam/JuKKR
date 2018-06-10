SUBROUTINE phicalc(iatom,lphi,visp,ipan,ircut,r,drdi,z,  &
        erefldau,phi,nspin,nsra)
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

IMPLICIT NONE
INCLUDE 'inc.p'
INTEGER NPOTD
PARAMETER (NPOTD= (2*KREL + (1-KREL)*NSPIND)*NATYPD)
DOUBLE PRECISION CVLIGHT
PARAMETER (CVLIGHT=274.0720442D0)

DOUBLE COMPLEX PHI(IRMD),PZ(IRMD,0:LMAXD),FZ(IRMD,0:LMAXD)
DOUBLE COMPLEX HAMF(IRMD,0:LMAXD),MASS(IRMD),DLOGDP(0:LMAXD)

DOUBLE PRECISION R(IRMD,NATYPD),VISP(IRMD,NPOTD),Z(NATYPD)
DOUBLE PRECISION DRDI(IRMD,NATYPD)
DOUBLE PRECISION RS(IRMD,0:LMAXD),S(0:LMAXD),DROR(IRMD),VPOT(IRMD)
DOUBLE PRECISION EREFLDAU

INTEGER IPAN(NATYPD),IRCUT(0:IPAND,NATYPD) 
INTEGER IATOM,IPOT1,NSRA,NSPIN,IRC1,IPAN1
!.. l-value for LDA+U
INTEGER LPHI

INTEGER IR,L1
DOUBLE PRECISION WINT(IRMD),WNORM,CUTOFF(IRMD)
DOUBLE COMPLEX CNORM,EZ,CZERO

czero = DCMPLX(0.d0,0.d0)
ipan1 = ipan(iatom)
irc1 = ircut(ipan1,iatom)

ipot1 = (iatom-1)*nspin + 1 ! points at the correct potential,
! i.e., 1,3,5,.. for NSPIN=2.

! -> set VPOT = [ V(UP) + V(DN) ]/2  for IATOM
!    only for the spherical part.

CALL dcopy(irmd,visp(1,ipot1),1,vpot(1),1)
IF ( nspin == 2 ) THEN
  CALL daxpy(irmd,1.d0,visp(1,ipot1+1),1,vpot(1),1)
  CALL dscal(irmd,0.5D0,vpot,1)
END IF

! -> Prepare and call REGSOL

DO l1 = 0,lmaxd
  IF ( nsra == 2 ) THEN
    s(l1) = DSQRT(DBLE(l1*l1+l1+1)-4.0D0*z(iatom)*z(iatom) /(cvlight*cvlight))
    IF ( z(iatom) == 0.0D0 ) s(l1) = DBLE(l1)
  ELSE
    s(l1) = DBLE(l1)
  END IF
  rs(1,l1) = 0.0D0
  DO ir = 2,irmd
    rs(ir,l1) = r(ir,iatom)**s(l1)
  END DO
END DO

DO ir = 2,irc1
  dror(ir) = drdi(ir,iatom)/r(ir,iatom)
END DO
ez = DCMPLX(erefldau,0.d0)

! --> this call of regsol within a non-lda+u calculation

CALL regsol(cvlight,ez,nsra,dlogdp,fz,hamf,mass,pz,dror,  &
    r(1,iatom),s,vpot,z(iatom),ipan(iatom),ircut(0,iatom),  &
    0,-1,0D0,cutoff,irmd,ipand,lmaxd)

! The result of regsol is (R*r)/r**l (in non-sra and similar in sra)


! --> set current angular momentum as LPHI

IF ( nsra == 2 ) THEN
  DO ir = 2,irc1
    pz(ir,lphi) = pz(ir,lphi)*rs(ir,lphi)
    fz(ir,lphi) = fz(ir,lphi)/cvlight*rs(ir,lphi)/cvlight
  END DO
ELSE
  DO ir = 2,irc1
    pz(ir,lphi) = pz(ir,lphi)*rs(ir,lphi)
    fz(ir,lphi) = czero
  END DO
END IF

! -> Copy result to PHI (only large component)

CALL zcopy(irmd,pz(1,lphi),1,phi(1),1)

! -> Normalise in sphere:
!    Note: If we have normalised in cell instead of sphere, the choice
!    M=0 would be arbitrary. One would have different normalisation for
!    each M, if the cell has no cubic symmetry.  Here we normalise in
!    sphere, to avoid this inconsistency.

DO ir = 1,irc1
  wint(ir) = dreal( DCONJG(phi(ir)) * phi(ir) )
END DO


! --> integrate, normalisation factor = WNORM

CALL simpk(wint,wnorm,ipan1,ircut(0,iatom),drdi(1,iatom))

! --> normalise PZ,FZ to unit probability in WS cell

cnorm = 1.d0/DSQRT(wnorm)
CALL zscal(irc1,cnorm,phi,1)

RETURN
END SUBROUTINE phicalc
