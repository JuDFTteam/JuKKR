module mod_rhocore

contains

subroutine rhocore(nsra, ispin, nspin, i1, drdi, r, visp, a, b, zat, ircut, &
  rhoc, ecore, ncore, lcore, cscl, vtrel, btrel, rmrel, drdirel, r2drdirel, &
  zrel, jwsrel, irshift, ecorerel, nkcore, kapcore)
  ! ===================================================================
  ! RELATIVISTIC TREATMENT OF CORE ELECTRONS   July/2002
  ! SEE ROUTINE <DRVCORE> FOR A SHORT DESCRIPTION OF THE VARIABLES
  use :: mod_datatypes, only: dp
  use global_variables
   use mod_corel
   use mod_drvcore
  implicit none
  real (kind=dp) :: a, b, zat
  integer :: jwsrel, zrel, irshift
  integer :: i1, ispin, ncore, nspin, nsra

  real (kind=dp) :: drdi(irmd), ecore(20*(krel+1)), r(irmd), rhoc(irmd, 2), &
    visp(irmd)
  integer :: ircut(0:ipand), lcore(20*(krel+1))
  ! ===================================================================
  ! ..
  ! .. Local Scalars ..
  ! ..
  real (kind=dp) :: ecorerel(krel*20+(1-krel), 2)
  integer :: nkcore(20), kapcore(20*2)
  real (kind=dp) :: cscl(krel*lmaxd+1)
  real (kind=dp) :: vtrel(irmd*krel+(1-krel))
  real (kind=dp) :: btrel(irmd*krel+(1-krel))
  real (kind=dp) :: drdirel(irmd*krel+(1-krel)), r2drdirel(irmd*krel+(1-krel)) &
    , rmrel(irmd*krel+(1-krel))
  ! .. External Subroutines ..
  ! --------------------------------------------------------------
  ! ipr=0 : do not write state dependent information
  real (kind=dp) :: qc, qc1, rmax
  integer :: ipr, nr
  save :: qc
  ! ipr=1 : write something
  ! ipr=2 : write all (for debugging)
  external :: corel, drvcore
  ! --------------------------------------------------------------

  ! =======================================================================
  ! non/scalar-relativistic OR relativistic

  ipr = 0

  if (ispin==1) qc = 0.0e0_dp
  nr = ircut(1)
  rmax = r(nr)


  ! =======================================================================
  if (krel==0) then
    ! =======================================================================

    call corel(nsra, ipr, i1, rhoc(1,ispin), visp, ecore, lcore, ncore, drdi, &
      zat, qc1, a, b, ispin, nspin, nr, rmax, irmd)
    ! non/scalar-relativistic OR relativistic
    if (ipr/=0) write (1337, fmt=100) i1
    qc = qc + qc1
    if (ispin==nspin) write (1337, fmt=110) zat, qc
    ! =======================================================================
  else
    ! *********************************************************************
    call drvcore(ipr, i1, lcore, ncore, cscl, vtrel, btrel, rmrel, a, b, &
      drdirel, r2drdirel, zrel, jwsrel, irshift, rhoc, ecorerel, nkcore, &
      kapcore, ecore, lmaxd, irmd)
  end if
  ! * For KREL = 1 (relativistic mode)                                  *
  ! *                                                                   *
  ! *  NPOTD = 2 * NATYPD                                               *
  return
100 format (1x, 5('*'), ' core-relaxation for ', i3, 'th cell', ' was done ', &
    5('*'))
110 format (4x, 'nuclear charge  ', f10.6, 9x, 'core charge =   ', f10.6)
end subroutine rhocore

end module mod_rhocore
