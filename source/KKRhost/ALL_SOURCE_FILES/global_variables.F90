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

  implicit none

  integer :: kBdG = 0         !! Switch between normal (noSOC or SOC) calculations and supermatrix calculations 
                              !! for Bogoliubov-de-Gennes formalism
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
  integer :: lmmaxso          !! (KREL+KORBIT+1)*(LMAX+1)**2
  integer :: lpotd
  integer :: nchebd
  integer :: maxmshd = 30         !! maximal number of different k-meshes
  logical :: linterface = .false. !! use 2D or 3D mode, if True a matching with semi-inifinite surfaces must be performed


end module global_variables


#ifdef CPP_MPI
!-------------------------------------------------------------------------------
!> Summary: MPI Briadcast of global variables
!> Author: Jonathan Chico
!> Category: KKRhost, communication, initialization
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!> MPI broadcast routine for global variables (i.e. array dimensions etc.)
!-------------------------------------------------------------------------------
subroutine bcast_global_variables()
  use :: mpi
  use :: global_variables ! use all parameters
  use :: mod_mympi, only: master
  implicit none

  integer :: n !! number of paramters that are broadcasted
  integer, allocatable :: blocklen1(:), etype1(:) !! blocklength of variuables in derived data type and list of MPI datatypes
  integer :: ierr !! error status
  integer :: mympitype1 !! derived data type for collective communication
  integer (kind=mpi_address_kind), allocatable :: disp1(:) !! MPI addresses
  integer (kind=mpi_address_kind) :: base !! base address of first entry


  n = 60
  allocate (blocklen1(n), etype1(n), disp1(n), stat=ierr)
  if (ierr/=0) stop 'error allocating arrays in bcast_global_variables'

  call mpi_get_address(n, disp1(1), ierr)
  call mpi_get_address(irid, disp1(2), ierr)
  call mpi_get_address(krel, disp1(3), ierr)
  call mpi_get_address(nfund, disp1(4), ierr)
  call mpi_get_address(ipand, disp1(5), ierr)
  call mpi_get_address(ngshd, disp1(6), ierr)
  call mpi_get_address(ncleb, disp1(7), ierr)
  call mpi_get_address(knoco, disp1(8), ierr)
  call mpi_get_address(iemxd, disp1(9), ierr)
  call mpi_get_address(irnsd, disp1(10), ierr)
  call mpi_get_address(nmaxd, disp1(11), ierr)
  call mpi_get_address(ishld, disp1(12), ierr)
  call mpi_get_address(naclsd, disp1(13), ierr)
  call mpi_get_address(nspotd, disp1(14), ierr)
  call mpi_get_address(ntperd, disp1(15), ierr)
  call mpi_get_address(ntrefd, disp1(16), ierr)
  call mpi_get_address(nsheld, disp1(17), ierr)
  call mpi_get_address(ncelld, disp1(18), ierr)
  call mpi_get_address(nspind, disp1(19), ierr)
  call mpi_get_address(knosph, disp1(20), ierr)
  call mpi_get_address(korbit, disp1(21), ierr)
  call mpi_get_address(kpoibz, disp1(22), ierr)
  call mpi_get_address(wlength, disp1(23), ierr)
  call mpi_get_address(nprincd, disp1(24), ierr)
  call mpi_get_address(nlayerd, disp1(25), ierr)
  call mpi_get_address(natomimpd, disp1(26), ierr)
  call mpi_get_address(lmaxd, disp1(27), ierr)
  call mpi_get_address(lmmaxd, disp1(28), ierr)
  call mpi_get_address(lmgf0d, disp1(29), ierr)
  call mpi_get_address(alm, disp1(30), ierr)
  call mpi_get_address(almgf0, disp1(31), ierr)
  call mpi_get_address(ndim_slabinv, disp1(32), ierr)
  call mpi_get_address(nembd, disp1(33), ierr)
  call mpi_get_address(nembd1, disp1(34), ierr)
  call mpi_get_address(nembd2, disp1(35), ierr)
  call mpi_get_address(nrd, disp1(36), ierr)
  call mpi_get_address(lm2d, disp1(37), ierr)
  call mpi_get_address(nclsd, disp1(38), ierr)
  call mpi_get_address(mmaxd, disp1(39), ierr)
  call mpi_get_address(npotd, disp1(40), ierr)
  call mpi_get_address(lmxspd, disp1(41), ierr)
  call mpi_get_address(lassld, disp1(42), ierr)
  call mpi_get_address(irmind, disp1(43), ierr)
  call mpi_get_address(nofgij, disp1(44), ierr)
  call mpi_get_address(nspindd, disp1(45), ierr)
  call mpi_get_address(nsatypd, disp1(46), ierr)
  call mpi_get_address(nrefd, disp1(47), ierr)
  call mpi_get_address(irmd, disp1(48), ierr)
  call mpi_get_address(naezd, disp1(49), ierr)
  call mpi_get_address(natypd, disp1(50), ierr)
  call mpi_get_address(lmpotd, disp1(51), ierr)
  call mpi_get_address(ntotd, disp1(52), ierr)
  call mpi_get_address(nrmaxd, disp1(53), ierr)
  call mpi_get_address(lmmaxso, disp1(54), ierr)
  call mpi_get_address(lpotd, disp1(55), ierr)
  call mpi_get_address(nchebd, disp1(56), ierr)
  call mpi_get_address(maxmshd, disp1(57), ierr)
  call mpi_get_address(kBdG, disp1(58), ierr)
  call mpi_get_address(linterface, disp1(59), ierr)
  call mpi_get_address(lnc, disp1(60), ierr)

  ! find displacements of variables
  base = disp1(1)
  disp1 = disp1 - base

  ! set length of variables in derived data type
  blocklen1(1:n) = 1

  ! set datatype of variables
  etype1(1:n-2) = mpi_integer
  etype1(n-1:n) = mpi_logical

  ! create new Type structure for derived data type
  call mpi_type_create_struct(n, blocklen1, disp1, etype1, mympitype1, ierr)
  if (ierr/=mpi_success) stop 'Problem in create_mpimask_t_inc'

  ! commit new type
  call mpi_type_commit(mympitype1, ierr)
  if (ierr/=mpi_success) stop 'error commiting create_mpimask_t_inc'

  ! broadcast derived data type

  call mpi_bcast(n, 1, mympitype1, master, mpi_comm_world, ierr)
  if (ierr/=mpi_success) stop 'error brodcasting t_inc'

  ! finally free auxiliary type and deallocate working arrays
  call mpi_type_free(mympitype1, ierr)
  deallocate (blocklen1, etype1, disp1, stat=ierr)
  if (ierr/=0) stop 'error deallocating arrays in bcast_global_variables'

end subroutine bcast_global_variables
#endif

