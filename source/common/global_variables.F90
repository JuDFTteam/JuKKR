!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------
!> Summary: Module containing the variables which were previously in the inc.p file
!> Author: Jonathan Chico
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!-------------------------------------------------------------------------------
module global_variables

  use mod_datatypes, only: dp
  implicit none

  integer :: kBdG = 0         !! Switch between normal (noSOC or SOC) calculations and supermatrix calculations for Bogoliubov-de-Gennes formalism
  integer :: irid = 350       !! Shape functions parameters in non-spherical part
  integer :: krel = 0         !! Switch for non-relativistic/relativistic (0/1) program. 
                              !! Attention: several other parameters depend explicitly on KREL,
                              !! they are set automatically Used for Dirac solver in ASA
  integer :: nfund = 289      !! Shape functions parameters in non-spherical part
  integer :: ipand = 50       !! Number of panels in non-spherical part
  integer :: ngshd = 60000    !! Shape functions parameters in non-spherical part
  integer :: ncleb = 784      !! Number of Clebsch-Gordon coefficients
  integer :: knoco = 0        !! (0/1) Collinear/Non-collinear magnetism (even in non-relativistic non-spin-orbit case)
  integer :: iemxd = 101      !! Dimension for energy-dependent arrays
  integer :: irnsd = 890      !! Number of radial mesh points in (RMT,...,RWS)
  integer :: nmaxd = 2000000  !! Paremeters for the Ewald summations
  integer :: ishld = 200000   !! Paremeters for the Ewald summations
  integer :: naclsd = 500     !! Maximum number of atoms in a TB-cluster
  integer :: nspotd = 2       !! Number of potentials for storing non-sph. potentials
  integer :: ntperd = 1       !! Parameter in broyden subroutines
  integer :: ntrefd = 0       !! parameter in broyden subroutine MUST BE 0 for the host program
  integer :: nsheld = 301     !! Number of blocks of the GF matrix that need to be calculated (NATYPD + off-diagonals in case of impurity)
  integer :: ncelld = 1       !! Number of cells (shapes) in non-spherical part
  integer :: nspind = 2       !! KREL+(1-KREL)*(NSPIN+1)
  integer :: knosph = 1       !! switch for spherical/non-spherical (0/1) program.
  integer :: korbit = 0       !! Spin-orbit/non-spin-orbit (1/0) added to the Schroedinger or SRA equations.
                              !! Works with FP. KREL and KORBIT cannot be both non-zero.
  integer :: kpoibz = 250000  !! Number of reciprocal space vectors

#ifdef __GFORTRAN__
  integer :: wlength = 4      !! Word length for direct access files, compiler dependent e.g. gfortran=4
#else
  integer :: wlength = 1      !! Word length for direct access files, compiler dependent e.g. ifort=1
#endif

  integer :: nprincd = 1      !! Number of principle layers, set to a number >= NRPINC in output of main0
  integer :: nlayerd = 1      !! Number of principal layers (NAEZD/NPRINCD) used in the inversion routines (independent on NATYPD)
  integer :: natomimpd = 150  !! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
  logical :: lnc = .true.     !! Coupled equations in two spins (switches true if KREL=1 or KORBIT=1 or KNOCO=1)

  integer :: lmaxd = 3        !! lmax cutoff
  integer :: lmmaxd = 16      !! (KREL+KORBIT+1)*(LMAXD+1)^2
  integer :: lmgf0d           !! (lmaxd+1)^2dimension of the reference system Green function, 
                              !! set up in the spin-independent non-relativstic (l,m_l)-representation
  integer :: alm              !! naez*lmmaxd
  integer :: almgf0           !! naezd*lmgf0d (<=alm)
  integer :: ndim_slabinv     !! nprincd*lmmaxd
  integer :: nembd = 20       !!
  integer :: nembd1 = 21      !! NEBMD+1
  integer :: nembd2 = 21      !! NEBMD+NAEZ
  integer :: nrd = 20000      !! Number of real space vectors
  integer :: lm2d
  integer :: nclsd
  integer :: mmaxd
  integer :: npotd
  integer :: lmxspd
  integer :: lassld
  integer :: irmind
  integer :: nofgij = 2       !! NATOMIMPD*NATOMIMPD+1 probably the same variable than NOFGIJD
  integer :: nspindd = 1      !! NSPIND-KORBIT
  integer :: nsatypd   = 1    !! (NATYPD-1)*KNOSPH+1
  integer :: nrefd = 1
  integer :: irmd = 900       !! Number of radial mesh points in (0,...,RWS)
  integer :: naezd = 1
  integer :: natypd = 1
  integer :: lmpotd
  integer :: ntotd = 80       !! IPAND+30
  integer :: nrmaxd
  integer :: lpotd
  integer :: nchebd
  integer :: maxmshd = 30         !! maximal number of different k-meshes
  logical :: linterface = .false. !! use 2D or 3D mode, if True a matching with semi-inifinite surfaces must be performed
  real (kind=dp) :: delta_BdG = 10**-4 !! initial value of Delta_BdG in Ry
  real (kind=dp) :: pot_ns_cutoff = -1 !! threshold below which non-spherical part of the potential is cut off (done in main2)


end module global_variables


