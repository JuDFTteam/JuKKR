!      MODULE MOD_RENORM_LLY
!      CONTAINS
      SUBROUTINE RENORMRHO_LLY(CDOS_LLY,RHOSPHER,IRMAX,IELAST,NSPIN, ! LLY Lloyd
     &         NATYP,DEN,LMAXP1,CONC,IESTART,IEEND,THETAS,NTCELL,RHO2NS)
! Renormalize the valence charge according to Lloyd's formula
      implicit none
      INCLUDE 'inc.p'
      INTEGER LMAXD1
      PARAMETER (LMAXD1=LMAXD+1)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER NPOTD
      PARAMETER (NPOTD= (2*(KREL+KORBIT) + 
     +           (1-(KREL+KORBIT))*NSPIND)*NATYPD)
! Input:
      INTEGER LMAXP1,NATYP,NSPIN,IRMAX         ! LMAX+1 for FP or LMAX for ASA, number of atoms,radial pts
      INTEGER IESTART,IEEND,IELAST             ! Starting and ending energy point for renormalization
      INTEGER NTCELL(NATYPD),NFU(NATYPD),LLMSP(NATYPD,NFUND) ! Shape-function details
      REAL*8 CONC(NATYPD)                      ! Concentration (for cpa)
      REAL*8 RHOSPHER(IRMAX,IELAST,NSPIN,NATYP)! Spherical component of normalized density/atom/energy/spin
      REAL*8 THETAS(IRID,NFUND,NCELLD)
      COMPLEX*16 DEN(0:LMAXD1,IEMXD,NPOTD)     ! Non-renormalized charge per atom
      COMPLEX*16 CDOS_LLY(IEMXD,NSPIND)        ! DOS according to Lloyd's formula
! Input/Output: 
      REAL*8 RHO2NS(IRMD,LMPOTD,NATYPD,2)      ! Valence charge density to be renormalized in its 
                                               ! spherical component
! Internal:
      INTEGER LL,IE,I1,ISPIN,IPOT,SPINDEGEN
      REAL*8 DENLOC(2),DENLLY(2)               ! Density from local summation and from Lloyd's formula
      REAL*8 RENORM(2)                         ! Renormalization constant for charge and spin density
      REAL*8 PI

      PI = 4.D0 * DATAN(1.D0)

      SPINDEGEN = 3 - NSPIN ! Spin degeneracy, 2 if nspin=1, 1 if nspin=2

      RENORM(:) = 1.D0
      DENLOC(:) = 0.D0
      DENLLY(:) = 0.D0


! Delete old spherical density
      RHO2NS(1:IRMD,1,1:NATYPD,1:2) = 0.D0

      DO IE = IESTART,IEEND
         DO ISPIN = 1,NSPIN
            DENLLY(ISPIN) = DIMAG(CDOS_LLY(IE,ISPIN)) * SPINDEGEN
            DENLOC(ISPIN) = 0.D0
            DO I1 = 1,NATYP
               IPOT = (I1-1) * NSPIN + ISPIN
               DO LL = 0,LMAXP1
                  DENLOC(ISPIN) = DENLOC(ISPIN) - 2.0D0 * CONC(I1) * 
     &                 DIMAG(DEN(LL,IE,IPOT))/PI/DBLE(NSPIN)
               END DO
            ENDDO
         ENDDO

         RENORM(1) = (DENLLY(1) + DENLLY(2)) / (DENLOC(1) + DENLOC(2)) ! Charge density
         RENORM(2) = (DENLLY(2) - DENLLY(1)) / (DENLOC(1) + DENLOC(2)) ! Spin density
         IF (LNC.OR.NSPIN.EQ.1) RENORM(2) = RENORM(1) ! Spins are coupled, only charge density is given by Lloyd's formula

         DO ISPIN = 1,NSPIN
            RHO2NS(1:IRMD,1,1:NATYP,ISPIN) = 
     &           RHO2NS(1:IRMD,1,1:NATYP,ISPIN) + 
     &           RHOSPHER(1:IRMD,IE,ISPIN,1:NATYP) * RENORM(ISPIN)  ! Integration weight included in RHOSPHER
         ENDDO

         WRITE(*,FMT='(A12,I5,2F16.12)') 
     &           'RENORM_LLY: ',IE,RENORM(1),RENORM(2)


      ENDDO


      END SUBROUTINE RENORM_LLY
!      END MODULE MOD_RENORM_LLY
