SUBROUTINE renormrho_lly(cdos_lly,rhospher,irmax,ielast,nspin, ! LLY Lloyd  &
    natyp,den,lmaxp1,conc,iestart,ieend,thetas,ntcell,rho2ns)
! Renormalize the valence charge according to Lloyd's formula
IMPLICIT NONE
INCLUDE 'inc.p'
INTEGER :: lmaxd1
PARAMETER (lmaxd1=lmaxd+1)
INTEGER :: lmpotd
PARAMETER (lmpotd= (lpotd+1)**2)
INTEGER :: npotd
PARAMETER (npotd= (2*(krel+korbit) + (1-(krel+korbit))*nspind)*natypd)
! Input:
INTEGER :: lmaxp1,natyp,nspin,irmax         ! LMAX+1 for FP or LMAX for ASA, number of atoms,radial pts
INTEGER :: iestart,ieend,ielast             ! Starting and ending energy point for renormalization
INTEGER :: ntcell(natypd),nfu(natypd),llmsp(natypd,nfund) ! Shape-function details
REAL*8 conc(natypd)                      ! Concentration (for cpa)
REAL*8 rhospher(irmax,ielast,nspin,natyp)! Spherical component of normalized density/atom/energy/spin
REAL*8 thetas(irid,nfund,ncelld)
COMPLEX*16 den(0:lmaxd1,iemxd,npotd)     ! Non-renormalized charge per atom
COMPLEX*16 cdos_lly(iemxd,nspind)        ! DOS according to Lloyd's formula
! Input/Output:
REAL*8 rho2ns(irmd,lmpotd,natypd,2)      ! Valence charge density to be renormalized in its spherical component
! Internal:
INTEGER :: ll,ie,i1,ispin,ipot,spindegen
REAL*8 denloc(2),denlly(2)               ! Density from local summation and from lloyd's formula
REAL*8 renorm(2)                         ! Renormalization constant for charge and spin density
REAL*8 pi

pi = 4.d0 * DATAN(1.d0)

spindegen = 3 - nspin ! Spin degeneracy, 2 if nspin=1, 1 if nspin=2

renorm(:) = 1.d0
denloc(:) = 0.d0
denlly(:) = 0.d0


! Delete old spherical density
rho2ns(1:irmd,1,1:natypd,1:2) = 0.d0

DO ie = iestart,ieend
  DO ispin = 1,nspin
    denlly(ispin) = DIMAG(cdos_lly(ie,ispin)) * spindegen
    denloc(ispin) = 0.d0
    DO i1 = 1,natyp
      ipot = (i1-1) * nspin + ispin
      DO ll = 0,lmaxp1
        denloc(ispin) = denloc(ispin) - 2.0D0 * conc(i1) *  &
            DIMAG(den(ll,ie,ipot))/pi/DBLE(nspin)
      END DO
    END DO
  END DO
  
  renorm(1) = (denlly(1) + denlly(2)) / (denloc(1) + denloc(2)) ! Charge density
  renorm(2) = (denlly(2) - denlly(1)) / (denloc(1) + denloc(2)) ! Spin density
  IF (lnc.OR.nspin == 1) renorm(2) = renorm(1) ! Spins are coupled, only charge density is given by lloyd's formula
  
  DO ispin = 1,nspin
    rho2ns(1:irmd,1,1:natyp,ispin) = rho2ns(1:irmd,1,1:natyp,ispin) +  &
        rhospher(1:irmd,ie,ispin,1:natyp) * renorm(ispin)  ! Integration weight included in rhospher
  END DO
  
  WRITE(1337,FMT='(A12,I5,2F16.12)') 'RENORM_LLY: ',ie,renorm(1),renorm(2)
  
  
END DO


END SUBROUTINE renorm_lly
!      END MODULE MOD_RENORM_LLY
