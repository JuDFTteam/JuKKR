!-------------------------------------------------------------------------------
! MODULE: global_variables
!> @brief Module containing the variables which were previously in the inc.p
!> file
!> @author Jonathan Chico
!-------------------------------------------------------------------------------
module global_variables

  implicit none

  integer :: irid !< Shape functions parameters in non-spherical part
  integer :: krel !< Switch for non-relativistic/relativistic (0/1) program. Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
  integer :: nfund !< Shape functions parameters in non-spherical part
  integer :: ipand !< Number of panels in non-spherical part
  integer :: ngshd !< Shape functions parameters in non-spherical part
  integer :: ncleb !< Number of Clebsch-Gordon coefficients
  integer :: knoco !< (0/1) Collinear/Non-collinear magnetism (even in non-relativistic non-spin-orbit case)
  integer :: iemxd !< Dimension for energy-dependent arrays
  integer :: irnsd
  integer :: nmaxd !< Paremeters for the Ewald summations
  integer :: ishld !< Paremeters for the Ewald summations
  integer :: naclsd !< Maximum number of atoms in a TB-cluster
  integer :: nspotd !< Number of potentials for storing non-sph. potentials
  integer :: ntperd !< Parameter in broyden subroutines
  integer :: ntrefd !< parameter in broyden subroutine MUST BE 0 for the host program
  integer :: nsheld !< Number of blocks of the GF matrix that need to be calculated (NATYPD + off-diagonals in case of impurity)
  integer :: ncelld !< Number of cells (shapes) in non-spherical part
  integer :: nspind !< KREL+(1-KREL)*(NSPIN+1)
  integer :: knosph !< switch for spherical/non-spherical (0/1) program.
  integer :: korbit !< Spin-orbit/non-spin-orbit (1/0) added to the Schroedinger or SRA equations. Works with FP. KREL and KORBIT cannot be both non-zero.
  integer :: kpoibz !< Number of reciprocal space vectors
  integer :: wlength !< Word length for direct access files, compiler dependent ifort/others (1/4)
  integer :: nprincd !< Number of principle layers, set to a number >= NRPINC in output of main0
  integer :: nlayerd !< Number of principal layers (NAEZD/NPRINCD) used in the inversion routines (independent on NATYPD)
  integer :: natomimpd !< Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
  logical :: lnc !< Coupled equations in two spins (switches true if KREL=1 or KORBIT=1 or KNOCO=1)

end module
