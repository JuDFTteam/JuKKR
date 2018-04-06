!-------------------------------------------------------------------------------
! MODULE: global_variables
!> @brief Module containing the variables which were previously in the inc.p
!> file
!> @author Jonathan Chico
!-------------------------------------------------------------------------------
module global_variables

   implicit none

   integer :: IRID      !< Shape functions parameters in non-spherical part
   integer :: KREL      !< Switch for non-relativistic/relativistic (0/1) program. Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
   integer :: NFUND     !< Shape functions parameters in non-spherical part
   integer :: IPAND     !< Number of panels in non-spherical part
   integer :: NGSHD     !< Shape functions parameters in non-spherical part
   integer :: NCLEB     !< Number of Clebsch-Gordon coefficients
   integer :: KNOCO     !< (0/1) Collinear/Non-collinear magnetism (even in non-relativistic non-spin-orbit case)
   integer :: IEMXD     !< Dimension for energy-dependent arrays
   integer :: IRNSD
   integer :: NMAXD     !< Paremeters for the Ewald summations
   integer :: ISHLD     !< Paremeters for the Ewald summations
   integer :: NACLSD    !< Maximum number of atoms in a TB-cluster
   integer :: NSPOTD    !< Number of potentials for storing non-sph. potentials
   integer :: NTPERD    !< Parameter in broyden subroutines
   integer :: NTREFD    !< parameter in broyden subroutine MUST BE 0 for the host program
   integer :: NSHELD    !< Number of blocks of the GF matrix that need to be calculated (NATYPD + off-diagonals in case of impurity)
   integer :: NCELLD    !< Number of cells (shapes) in non-spherical part
   integer :: NSPIND    !< KREL+(1-KREL)*(NSPIN+1)
   integer :: KNOSPH    !< switch for spherical/non-spherical (0/1) program.
   integer :: KORBIT    !< Spin-orbit/non-spin-orbit (1/0) added to the Schroedinger or SRA equations. Works with FP. KREL and KORBIT cannot be both non-zero.
   integer :: KPOIBZ    !< Number of reciprocal space vectors
   integer :: WLENGTH   !< Word length for direct access files, compiler dependent ifort/others (1/4)
   integer :: NPRINCD   !< Number of principle layers, set to a number >= NRPINC in output of main0
   integer :: NLAYERD   !< Number of principal layers (NAEZD/NPRINCD) used in the inversion routines (independent on NATYPD)
   integer :: NATOMIMPD !< Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
   logical :: LNC       !< Coupled equations in two spins (switches true if KREL=1 or KORBIT=1 or KNOCO=1)

end module global_variables
