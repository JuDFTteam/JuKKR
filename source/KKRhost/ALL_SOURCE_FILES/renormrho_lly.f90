subroutine renormrho_lly(cdos_lly, rhospher, irmax, ielast, nspin, ! LLY Lloyd  &

  natyp, den, lmaxp1, conc, iestart, ieend, thetas, ntcell, rho2ns)
! Renormalize the valence charge according to Lloyd's formula
      Use mod_datatypes, Only: dp
  implicit none
  include 'inc.p'
  integer :: lmaxd1
  parameter (lmaxd1=lmaxd+1)
  integer :: lmpotd
  parameter (lmpotd=(lpotd+1)**2)
  integer :: npotd
  parameter (npotd=(2*(krel+korbit)+(1-(krel+korbit))*nspind)*natypd)
! Concentration (for cpa)
  integer :: lmaxp1, natyp, nspin, irmax ! Spherical component of normalized density/atom/energy/spin
  integer :: iestart, ieend, ielast ! Non-renormalized charge per atom
  integer :: ntcell(natypd), nfu(natypd), llmsp(natypd, nfund) ! DOS according to Lloyd's formula
  real (kind=dp) :: conc(natypd) ! Input/Output:
  real (kind=dp) :: rhospher(irmax, ielast, nspin, natyp) ! Valence charge density to be renormalized in its spherical component
  real (kind=dp) :: thetas(irid, nfund, ncelld)
  complex (kind=dp) :: den(0:lmaxd1, iemxd, npotd) ! Internal:
  complex (kind=dp) :: cdos_lly(iemxd, nspind) ! Density from local summation and from lloyd's formula
! Renormalization constant for charge and spin density
  real (kind=dp) :: rho2ns(irmd, lmpotd, natypd, 2) 

  integer :: ll, ie, i1, ispin, ipot, spindegen
  real (kind=dp) :: denloc(2), denlly(2) ! Spin degeneracy, 2 if nspin=1, 1 if nspin=2
  real (kind=dp) :: renorm(2) 
  real (kind=dp) :: pi

  pi = 4.d0*datan(1.d0)

  spindegen = 3 - nspin ! Delete old spherical density

  renorm(:) = 1.d0
  denloc(:) = 0.d0
  denlly(:) = 0.d0

! Charge density
! Spin density
  rho2ns(1:irmd, 1, 1:natypd, 1:2) = 0.d0
! Spins are coupled, only charge density is given by lloyd's formula
  do ie = iestart, ieend
    do ispin = 1, nspin
      denlly(ispin) = dimag(cdos_lly(ie,ispin))*spindegen
      denloc(ispin) = 0.d0
      do i1 = 1, natyp
        ipot = (i1-1)*nspin + ispin
        do ll = 0, lmaxp1
          denloc(ispin) = denloc(ispin) - 2.0d0*conc(i1)*dimag(den(ll,ie,ipot) &
            )/pi/dble(nspin)
        end do
      end do
    end do

    renorm(1) = (denlly(1)+denlly(2))/(denloc(1)+denloc(2)) ! Integration weight included in rhospher
    renorm(2) = (denlly(2)-denlly(1))/(denloc(1)+denloc(2)) 
    if (lnc .or. nspin==1) renorm(2) = renorm(1) 

    do ispin = 1, nspin
      rho2ns(1:irmd, 1, 1:natyp, ispin) = rho2ns(1:irmd, 1, 1:natyp, ispin) + &
        rhospher(1:irmd, ie, ispin, 1:natyp)*renorm(ispin) 
    end do

    write (1337, fmt='(A12,I5,2F16.12)') 'RENORM_LLY: ', ie, renorm(1), &
      renorm(2)
! LLY Lloyd  &
! Renormalize the valence charge according to Lloyd's formula
  end do
! set to 1 if NEWSOSOL under RUNOPT, otherwise 0
! SET ACCORDING TO lmax VALUE OF INPUTCARD
end subroutine
!      PARAMETER ( NRD = 20000, KPOIBZ = 32000 )
