  subroutine init_param()
! These values should be specified somewhere in the main program
! or read in from an input file
  use global

  implicit none

! Global parameters
  lmmax  = (nlmax+1)**2                  ! lm components for GFs
  nlmax2 = 2*nlmax
  lmmax2 = (2*nlmax+1)**2                ! L1 and L2 in Gaunt numbers
  nlmax4 = 4*nlmax
  lmmax4 = (4*nlmax+1)**2                ! L3 in Gaunt numbers
  nsmax2 = nsmax**2                      ! spin components; set to 2 if ANYTHING is magnetic
  nlms   = lmmax*nsmax                   ! total size of scattering matrices with spin
  nalms  = nasusc*nlms                   ! dimensions for structural GF
!***************************************************************************************************
! Allocate storage for input parameters and control information
  allocate(iasusc(nasusc))                      ! which atoms for projection
  iasusc = -1
  allocate(iasusc2(nasusc))                     ! which atoms for susceptibility
  iasusc2 = -1
  allocate(igroup(nasusc))                      ! how many atoms in each group
  iasusc2 = -1
! --------------------------------------------------------------------
  allocate(inobxc(nasusc))                      ! atom-dependent spin density killer
  inobxc = 0
! --------------------------------------------------------------------
  allocate(iadmat(nasusc))                      ! atom-dependent density matrix switches
  allocate(ildmat(0:nlmax,nsmax,nasusc))        ! atom-dependent density matrix switches
  iadmat = 0; ildmat = 0
! --------------------------------------------------------------------
  allocate(isoc(nasusc))                        ! which atoms for SOC
  allocate(socscaling(nasusc))                  ! scaling of SOC for each atom
  isoc = 0; socscaling = 0.d0
! --------------------------------------------------------------------
  allocate(ildau(nasusc))                       ! which atoms for LDA+U
  allocate(ueff(0:nlmax,nasusc))                ! Dudarev U parameters
  allocate(jeff(0:nlmax,nasusc))                ! Dudarev J parameters
  ildau = 0; ueff = 0.d0; jeff = 0.d0
! --------------------------------------------------------------------
  allocate(ijij(nasusc))                        ! which atoms for Jij
  ijij = 0
! --------------------------------------------------------------------
  allocate(ibfield(nasusc))                     ! which atoms for B field
  allocate(blen(nasusc),bdir(3,nasusc))         ! strength and direction of B field for each atom
  allocate(bconlen(nasusc),bcondir(3,nasusc))   ! strength and direction of constraining B field for each atom
  ibfield = 0; blen = 0.d0; bdir = 0.d0
  bconlen = 0.d0; bcondir = 0.d0
! --------------------------------------------------------------------
  allocate(isusc(nasusc))                       ! susc options
  allocate(ikha(nasusc))                        ! Hartree kernel
  allocate(ikxc(nasusc))                        ! xc options
  isusc = 0; ikha = 0; ikxc = 0
! --------------------------------------------------------------------
  allocate(iwsusc(0:nlmax,nsmax,nasusc))        ! which l-blocks to use for GF
  allocate(issusc(nasusc))                      ! atom magnetic or not
  allocate(ewsusc(nbmax,0:nlmax,nsmax,nasusc))  ! projection energies -- nbmax should be fine from itc = 1
  allocate(nowfns(0:nlmax,nsmax,nasusc))        ! were all wfns stored?
  allocate(nobasis(0:nlmax,nsmax,nasusc))       ! was a basis constructed?
  iwsusc = -1; issusc = 0; ewsusc = 0.d0
  nowfns = .true.; nobasis = .true.
! ----------------------------------------------------------------------
! Pointers for spin labels: is=1 is down and is=2 is up
! Spin flip and then spin diagonals
  is2i(1,1) = 4; is2i(1,2) = 2
  is2i(2,1) = 1; is2i(2,2) = 3
  i2is(:,1) = (/2,1/); i2is(:,2) = (/1,2/)
  i2is(:,3) = (/2,2/); i2is(:,4) = (/1,1/)
! ----------------------------------------------------------------------
! Number of panel > 1
  allocate(npanat(nasusc))
  allocate(ircutat(npanmax+1,nasusc))
  npanat  = 0
  ircutat = 0

! Flag these things as initialized
  noparameters = .false.
! All done!
  end subroutine init_param
