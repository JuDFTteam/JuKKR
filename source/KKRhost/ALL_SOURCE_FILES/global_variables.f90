!-------------------------------------------------------------------------------
! MODULE: global_variables
!> @brief Module containing the variables which were previously in the inc.p
!> file
!> @author Jonathan Chico
!-------------------------------------------------------------------------------
    Module global_variables

      Implicit None

      Integer :: irid !< Shape functions parameters in non-spherical part
      Integer :: krel !< Switch for non-relativistic/relativistic (0/1) program. Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
      Integer :: nfund !< Shape functions parameters in non-spherical part
      Integer :: ipand !< Number of panels in non-spherical part
      Integer :: ngshd !< Shape functions parameters in non-spherical part
      Integer :: ncleb !< Number of Clebsch-Gordon coefficients
      Integer :: knoco !< (0/1) Collinear/Non-collinear magnetism (even in non-relativistic non-spin-orbit case)
      Integer :: iemxd !< Dimension for energy-dependent arrays
      Integer :: irnsd
      Integer :: nmaxd !< Paremeters for the Ewald summations
      Integer :: ishld !< Paremeters for the Ewald summations
      Integer :: naclsd !< Maximum number of atoms in a TB-cluster
      Integer :: nspotd !< Number of potentials for storing non-sph. potentials
      Integer :: ntperd !< Parameter in broyden subroutines
      Integer :: ntrefd !< parameter in broyden subroutine MUST BE 0 for the host program
      Integer :: nsheld !< Number of blocks of the GF matrix that need to be calculated (NATYPD + off-diagonals in case of impurity)
      Integer :: ncelld !< Number of cells (shapes) in non-spherical part
      Integer :: nspind !< KREL+(1-KREL)*(NSPIN+1)
      Integer :: knosph !< switch for spherical/non-spherical (0/1) program.
      Integer :: korbit !< Spin-orbit/non-spin-orbit (1/0) added to the Schroedinger or SRA equations. Works with FP. KREL and KORBIT cannot be both non-zero.
      Integer :: kpoibz !< Number of reciprocal space vectors
      Integer :: wlength !< Word length for direct access files, compiler dependent ifort/others (1/4)
      Integer :: nprincd !< Number of principle layers, set to a number >= NRPINC in output of main0
      Integer :: nlayerd !< Number of principal layers (NAEZD/NPRINCD) used in the inversion routines (independent on NATYPD)
      Integer :: natomimpd !< Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
      Logical :: lnc !< Coupled equations in two spins (switches true if KREL=1 or KORBIT=1 or KNOCO=1)

      integer :: lmaxd !< lmax cutoff
      integer :: lmmaxd !< (KREL+KORBIT+1)*(LMAXD+1)^2
      integer :: lmgf0d !< (lmaxd+1)^2
      integer :: alm !< naez*lmmaxd
      integer :: almgf0 !< naezd*lmgf0d (<=alm)
      integer :: ndim_slabinv !< nprincd*lmmaxd
      integer :: nembd1 !< NEBMD+1
      integer :: nembd2 !< NEBMD+NAEZ
      integer :: nrd
      integer :: lm2d
      integer :: nclsd
      integer :: mmaxd
      integer :: npotd
      integer :: lmxspd
      integer :: lassld
      integer :: irmind
      integer :: nofgij
      integer :: nspindd
      integer :: nsatypd
      integer :: nrefd
      integer :: irmd
      integer :: naezd
      integer :: natypd
      integer :: lmpotd
      integer :: ntotd
      integer :: nrmaxd
      integer :: nembd

      logical :: linterface

    End Module
