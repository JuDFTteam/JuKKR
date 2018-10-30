      SUBROUTINE PHICALC(IATOM,LPHI,VISP,IPAN,IRCUT,R,DRDI,Z,
     &                   EREFLDAU,PHI,NSPIN,NSRA)
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
      IMPLICIT NONE
      INCLUDE 'inc.p'
C     ..
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
C     .. l-value for LDA+U
      INTEGER LPHI
C
      INTEGER IR,L1
      DOUBLE PRECISION WINT(IRMD),WNORM,CUTOFF(IRMD)
      DOUBLE COMPLEX CNORM,EZ,CZERO
C
      CZERO = DCMPLX(0.D0,0.D0)
      IPAN1 = IPAN(IATOM)
      IRC1 = IRCUT(IPAN1,IATOM)
C
      IPOT1 = (IATOM-1)*NSPIN + 1 ! points at the correct potential,
                                  ! i.e., 1,3,5,.. for NSPIN=2.
C
C -> set VPOT = [ V(UP) + V(DN) ]/2  for IATOM
C    only for the spherical part.
C
      CALL DCOPY(IRMD,VISP(1,IPOT1),1,VPOT(1),1)
      IF ( NSPIN.EQ.2 ) THEN
          CALL DAXPY(IRMD,1.D0,VISP(1,IPOT1+1),1,VPOT(1),1)
          CALL DSCAL(IRMD,0.5D0,VPOT,1)
      END IF
C
C -> Prepare and call REGSOL
C
      DO L1 = 0,LMAXD
          IF ( NSRA.EQ.2 ) THEN
              S(L1) = DSQRT(DBLE(L1*L1+L1+1)-4.0D0*Z(IATOM)*Z(IATOM)
     &             /(CVLIGHT*CVLIGHT))
              IF ( Z(IATOM).EQ.0.0D0 ) S(L1) = DBLE(L1)
          ELSE
              S(L1) = DBLE(L1)
          END IF
          RS(1,L1) = 0.0D0
          DO IR = 2,IRMD
              RS(IR,L1) = R(IR,IATOM)**S(L1)
          END DO
      END DO
C     
      DO IR = 2,IRC1
          DROR(IR) = DRDI(IR,IATOM)/R(IR,IATOM)
      END DO
      EZ = DCMPLX(EREFLDAU,0.D0)
C
C --> this call of regsol within a non-lda+u calculation
C
      CALL REGSOL(CVLIGHT,EZ,NSRA,DLOGDP,FZ,HAMF,MASS,PZ,DROR,
     &            R(1,IATOM),S,VPOT,Z(IATOM),IPAN(IATOM),IRCUT(0,IATOM),
     &            0,-1,0D0,CUTOFF,IRMD,IPAND,LMAXD)
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
C --> integrate, normalisation factor = WNORM
C
      CALL SIMPK(WINT,WNORM,IPAN1,IRCUT(0,IATOM),DRDI(1,IATOM))
C
C --> normalise PZ,FZ to unit probability in WS cell
C
      CNORM = 1.D0/DSQRT(WNORM)
      CALL ZSCAL(IRC1,CNORM,PHI,1)
C      
      RETURN
      END
