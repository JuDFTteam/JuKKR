      MODULE MOD_PHICALC
      CONTAINS

      SUBROUTINE PHICALC(LPHI,PHI,VISP
     &     ,IPAN,IRCUT,R,DRDI,Z,EREFLDAU,IDOLDAU,WLDAUAV,CUTOFF
     &     ,IATOM,NSPIN,NSRA,LMAXD,IRMD)
C **********************************************************
C  Calculates test functions PHI for LDA+U.
C  Only spherical part of the potential is used.
C  PHI are then normalised within Voronoi Circumscribing Sphere.
C  In sra treatment, only large component is considered as 
C  important and normalised although small component is also 
C  calculated.
C
C  PHI is here one-dim. array (PHI(IRMD)), so the subroutine 
C  must be called for each atom seperately.
C **********************************************************
      USE mod_regsol
      USE mod_simpk
      implicit none
      REAL*8 CVLIGHT
      PARAMETER (CVLIGHT=274.0720442D0)

! Input
      INTEGER IATOM,NSRA,NSPIN,IRMD,LMAXD
      INTEGER LPHI                 ! l-value for LDA+U
      INTEGER IPAN,IRCUT(0:IPAN)   ! IPAN,IRCUT(0:IPAND)
      INTEGER IDOLDAU
      REAL*8 DRDI(:)             ! DRDI(IRMD,NATYPD)
      REAL*8 R(:),VISP(:,:),Z ! R(IRMD),VISP(IRMD,NPOTD)
      REAL*8 EREFLDAU,WLDAUAV(2),CUTOFF(:)

! Output:
      COMPLEX*16 PHI(:)          ! PHI(IRMD)

! Inside
      REAL*8,ALLOCATABLE :: RS(:,:),S(:),DROR(:)  ! RS(IRMD,0:LMAXD),S(0:LMAXD),DROR(IRMD)
      DOUBLE COMPLEX,ALLOCATABLE :: HAMF(:,:),MASS(:),DLOGDP(:)     ! HAMF(IRMD,0:LMAXD),MASS(IRMD),DLOGDP(0:LMAXD)
      REAL*8,ALLOCATABLE  :: VAVRG(:),WINT(:)    ! VAVRG(IRMD),WINT(IRMD)
      COMPLEX*16,ALLOCATABLE :: PZ(:,:),FZ(:,:) ! PZ(IRMD,0:LMAXD),FZ(IRMD,0:LMAXD)


      INTEGER IPOT1,IRS1,IRC1,IPAN1

      INTEGER ILM2,LM1,LM2,LM3,II,IR,L1,MMAX,LMAX
      REAL*8 WNORM,WLDAUAVUD
      COMPLEX*16 CNORM,EZ,CZERO



      ALLOCATE ( RS(IRMD,0:LMAXD),S(0:LMAXD),DROR(IRMD) )
      ALLOCATE( HAMF(IRMD,0:LMAXD),  MASS(IRMD),  DLOGDP(0:LMAXD) )
      ALLOCATE( VAVRG(IRMD), WINT(IRMD) )
      ALLOCATE( PZ(IRMD,0:LMAXD),FZ(IRMD,0:LMAXD) )

      CZERO = DCMPLX(0.D0,0.D0)
      MMAX = 2*LPHI + 1
      IRS1 = IRCUT(1)
      IRC1 = IRCUT(IPAN)

      IPOT1 = (IATOM-1)*NSPIN + 1 ! points at the correct potential,
                                  ! i.e., 1,3,5,.. for NSPIN=2.
C
C --> set VAVRG = [ V(UP) + V(DN) ]/2  for IATOM
C     only for the spherical part.
C
      WLDAUAVUD = WLDAUAV(1)
      CALL DCOPY(IRMD,VISP(1,IPOT1),1,VAVRG(1),1)
      IF ( NSPIN.EQ.2 ) THEN
          CALL DAXPY(IRMD,1.D0,VISP(1,IPOT1+1),1,VAVRG(1),1)
          CALL DSCAL(IRMD,0.5D0,VAVRG,1)
          WLDAUAVUD = 0.5D0 * (WLDAUAV(1) + WLDAUAV(2))
      END IF


C Prepare and call REGSOL
C
      DO L1 = 0,LMAXD
          IF ( NSRA.GE.2 ) THEN
              S(L1) = DSQRT(DBLE(L1*L1+L1+1)-4.0D0*Z*Z
     &             /(CVLIGHT*CVLIGHT))
              IF ( Z.EQ.0.0D0 ) S(L1) = DBLE(L1)
          ELSE
              S(L1) = DBLE(L1)
          END IF
          RS(1,L1) = 0.0D0
          DO IR = 2,IRC1
              RS(IR,L1) = R(IR)**S(L1)
          END DO
      END DO
C     
      DO IR = 2,IRC1
          DROR(IR) = DRDI(IR)/R(IR)
      END DO
      EZ = DCMPLX(EREFLDAU,0.D0)
      CALL REGSOL(CVLIGHT,EZ,NSRA,DLOGDP,FZ,HAMF,MASS,PZ,DROR,R,
     &     S,VAVRG,Z,IPAN,IRCUT,IDOLDAU,LPHI,WLDAUAVUD,CUTOFF,
     +     IRMD,IPAN,LMAXD)

C The result of regsol is (R*r)/r**l (in non-sra and similar in sra)
C
C
C --> set current angular momentum as LPHI
C
      IF ( NSRA.GE.2 ) THEN
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
C Copy result to PHI (only large component)
      CALL ZCOPY(IRMD,PZ(1,LPHI),1,PHI(1),1)

C Prepare Gaunt coefficients in array GSH and pass to array GAUNT:
c     CALL RINIT(LMMAXD*LMMAXD*LMPOTD,GAUNT(1,1,1))
c     CALL GAUNT2(WG,YRG)
c     CALL SHAPE(LMAXD,1,GSH,ILM,IMAXSH,LMSP,NTCELL,WG,YRG)
c     DO II = 1,IMAXSH(LMMAXD)
c         LM1 = ILM(II,1)
c         LM2 = ILM(II,2)
c         LM3 = ILM(II,3)
c         GAUNT(LM1,LM2,LM3) = GSH(II)
c     END DO
C
C --> prepare integrand for normalisation of the reg. wavefunctions
C     
C     WINT(r)= R(r)**2 * W2(r), W2(r)=Sum_{LLL'} C_{LLL'} Theta_{L'}(r)
c     CALL RINIT(IRMD,WINT(1))
c     CALL RINIT(IRMD,W2(1))

c     LM1 = LPHI**2+LPHI+1 !Set LM1 to (L=LPHI,M=0) (See Note below)
c     DO LM2 = 1,LMPOTD !Loop over shapes
c         IF (LMSP(IATOM,LM2).GT.0) THEN
c         DO IR = IRS1+1,IRC1
c             W2(IR) = W2(IR) + GAUNT(LM1,LM1,LM2) 
c    &                 * THETAS(IR-IRS1,LM2,ICELL)
c         ENDDO
c         END IF
c     ENDDO
C Normalise in cell:
c     DO IR = 1,IRS1
c         WINT(IR) = DREAL( DCONJG(PHI(IR)) * PHI(IR) )
c     END DO
c     DO IR = IRS1+1,IRC1
c         WINT(IR) = DREAL( DCONJG(PHI(IR)) * PHI(IR) )
c    &               * W2(IR)
c     ENDDO
C Or, Normalise in sphere:
      PHI(:) = CUTOFF(:) * PHI(:)
      DO IR = 1,IRC1
          WINT(IR) = DREAL( DCONJG(PHI(IR)) * PHI(IR) )
      ENDDO

C
C --> integrate, normalisation factor = WNORM
C
      CALL SIMPK(WINT,WNORM,IPAN,IRCUT,DRDI,IPAN)
C
C --> normalise PZ,FZ to unit probability in WS cell
C

      CNORM = 1.D0/DSQRT(WNORM)
      CALL ZSCAL(IRMD,CNORM,PHI,1)
   
      DEALLOCATE( RS,S,DROR )
      DEALLOCATE( HAMF,  MASS,  DLOGDP )
      DEALLOCATE( VAVRG, WINT )
      DEALLOCATE( PZ, FZ )
   
      
      RETURN
      END SUBROUTINE PHICALC

C     Note: If we would normalise in cell instead of sphere, the choice
C     M=0 would be arbitrary. One would have different normalisation for
C     each M, if the cell has no cubic symmetry.  Here we normalise in
C     sphere, to avoid this inconsistency.


      END MODULE MOD_PHICALC
