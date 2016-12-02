      SUBROUTINE LDAUPHI(LPHI,VISP,IPAN,IRCUT,R,DRDI,ZAT,
     &                   EREFLDAU,PHI,NSPIN,NSRA,
     &                   NLDAU,LLDAU,
C                        new input parameters after inc.p removal
     &                   lmaxd, irmd, ipand)
C *********************************************************************
C *                                                                   *
C * Calculates test functions PHI for LDA+U.                          *
C * Only spherical part of the potential is used.                     *
C * PHI are then normalised within Voronoi Circumscribing Sphere.     *
C * In SRA treatment, only large component is considered as important *
C * and normalised although small component is also calculated.       *
C *                                                                   *
C * PHI is here one-dime array -- PHI(IRMD) -- so the subroutine must *
C * be called for each atom seperately.                               *
C *                                  ph. mavropoulos, juelich 2004    *
C *                                                                   *
C *********************************************************************
      use SingleSiteHelpers_mod, only: regsol
      use Quadrature_mod, only: simpson
      IMPLICIT NONE

      INTEGER lmaxd
      INTEGER irmd
      INTEGER ipand

      DOUBLE PRECISION CVLIGHT
      PARAMETER        (CVLIGHT=274.0720442D0)
C
C global arrays ..
      DOUBLE COMPLEX   PHI(IRMD)
      DOUBLE PRECISION R(IRMD),
     +                 VISP(IRMD,2),
     +                 DRDI(IRMD)
      INTEGER          LLDAU(LMAXD + 1),IRCUT(0:IPAND)
C global scalars ..
      DOUBLE PRECISION EREFLDAU,ZAT
      INTEGER          IPAN,NSPIN,NSRA,LPHI,NLDAU
C ..
C local arrays ..
      DOUBLE COMPLEX   PZ(IRMD,0:LMAXD),FZ(IRMD,0:LMAXD),
     +                 HAMF(IRMD,0:LMAXD),MASS(IRMD),DLOGDP(0:LMAXD)
      DOUBLE PRECISION RS(IRMD,0:LMAXD),S(0:LMAXD),
     +                 DROR(IRMD),VPOT(IRMD),
     +                 WINT(IRMD),WMLDAUAV(LMAXD + 1),
     +                 LDAUCUT(IRMD)
C ..
C local scalars ..
      DOUBLE COMPLEX   CNORM,EZ,CZERO
      DOUBLE PRECISION WNORM
      INTEGER          IRC1,IR,L1,ILDAU
      LOGICAL          LDAUTMP
C ..
C
      CZERO    = DCMPLX(0.D0,0.D0)
      IRC1     = IRCUT(IPAN)
      DO ILDAU=1,NLDAU
        WMLDAUAV(ILDAU) = 0.0D0
      ENDDO
      DO IR=1,IRMD
        LDAUCUT(IR) = 1.0D0
      ENDDO
C
C
C -> set VPOT = [ V(UP) + V(DN) ]/2  for IATOM
C    only for the spherical part.
C
C
      CALL DCOPY(IRMD,VISP(1,1),1,VPOT(1),1)
      IF ( NSPIN.EQ.2 ) THEN
        CALL DAXPY(IRMD,1.D0,VISP(1,2),1,VPOT(1),1)
        CALL DSCAL(IRMD,0.5D0,VPOT,1)
      ENDIF
C
C -> Prepare and call REGSOL
C
      DO L1 = 0,LMAXD
          IF ( NSRA.EQ.2 ) THEN
              S(L1) = DSQRT(DBLE(L1*L1+L1+1)-4.0D0*ZAT*ZAT
     &             /(CVLIGHT*CVLIGHT))
              IF ( ZAT.EQ.0.0D0 ) S(L1) = DBLE(L1)
          ELSE
              S(L1) = DBLE(L1)
          END IF
          RS(1,L1) = 0.0D0
          DO IR = 2,IRMD
              RS(IR,L1) = R(IR)**S(L1)
          END DO
      END DO
C     
      DO IR = 2,IRC1
          DROR(IR) = DRDI(IR)/R(IR)
      END DO
      EZ = DCMPLX(EREFLDAU,0.D0)
C
C --> this call of regsol within a non-lda+u calculation
C
!       CALL REGSOL(CVLIGHT,EZ,NSRA,DLOGDP,FZ,HAMF,MASS,PZ,DROR,
!      &            R(1),S,VPOT,ZAT,IPAN,IRCUT(0),
!      &            0,-1,0D0,CUTOFF,IRMD,IPAND,LMAXD)
C
      LDAUTMP = .FALSE.
C
      CALL REGSOL(CVLIGHT,EZ,NSRA,DLOGDP,FZ,HAMF,MASS,PZ,DROR,
     &            R(1),S,VPOT,ZAT,IPAN,IRCUT(0),
     &            IRMD,IPAND,LMAXD,
     &            LDAUTMP,NLDAU,LLDAU,WMLDAUAV,LDAUCUT)
C
C The result of regsol is (R*r)/r**l (in non-sra and similar in sra)
C
C
C --> set current angular momentum as LPHI
C
      IF ( NSRA.EQ.2 ) THEN
          DO IR = 2,IRC1
              PZ(IR,LPHI) = PZ(IR,LPHI)*RS(IR,LPHI)
              FZ(IR,LPHI) = FZ(IR,LPHI)/CVLIGHT*RS(IR,LPHI)/CVLIGHT
          END DO
      ELSE
          DO IR = 2,IRC1
              PZ(IR,LPHI) = PZ(IR,LPHI)*RS(IR,LPHI)
              FZ(IR,LPHI) = CZERO
          END DO
      END IF
C
C -> Copy result to PHI (only large component)
C
      CALL ZCOPY(IRMD,PZ(1,LPHI),1,PHI(1),1)
C
C
C -> Normalise in sphere:
C    Note: If we have normalised in cell instead of sphere, the choice
C    M=0 would be arbitrary. One would have different normalisation for
C    each M, if the cell has no cubic symmetry.  Here we normalise in
C    sphere, to avoid this inconsistency.
C
      DO IR = 1,IRC1
          WINT(IR) = DREAL( DCONJG(PHI(IR)) * PHI(IR) )
      ENDDO
C
C
C --> integrate, normalisation factor = WNORM
C
c     CALL SIMPK(WINT,WNORM,IPAN,IRCUT(0),DRDI(1))
      WNORM = simpson(WINT, IPAN, IRCUT(0:), DRDI(1:))
C
C --> normalise PZ,FZ to unit probability in WS cell
C
      CNORM = 1.D0/DSQRT(WNORM)
      CALL ZSCAL(IRC1,CNORM,PHI,1)
C
C
      RETURN
      END
